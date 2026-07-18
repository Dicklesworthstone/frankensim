//! fs-ledger: the Design Ledger v0 (plan §11.2, Bet 10; Decalogue P4/P9).
//!
//! FrankenSQLite-backed system of record: content-addressed artifacts,
//! event-sourced ops with the frozen Five Explicits, lineage edges, metric
//! time series, the autotuner cache, and the fine-grained event stream —
//! one file per project, WAL mode, snapshot-isolated readers.
//!
//! Layer: L6 (HELM). Runtime dependencies: `std` + `fsqlite` (constellation).
//!
//! Concurrency contract (from FrankenSQLite's documented model): open one
//! [`Ledger`] (one connection) per thread within one process. Readers get
//! snapshot isolation and never block the appending writer; contention
//! surfaces as a retryable [`LedgerError::Busy`], never a hang or a silent
//! lost write. Multi-process multi-writer use is a no-claim boundary
//! (CONTRACT.md).
//!
//! Time travel, forkable worlds, `explain()`, the replay audit, and GC
//! live in the [`travel`] module (schema v2).

pub mod colors;
pub mod crosswalk;
pub mod hash;
pub mod identity_migration;
pub mod schema;
pub mod session_registry;
pub mod state_checkpoint;
pub mod tombstone;
pub mod travel;
pub mod vcs;

pub use colors::{
    COLOR_DEMOTION_ROW_SCHEMA_VERSION, COLOR_WRITE_ROW_SCHEMA_VERSION, ColorDemotion, ColorGraph,
    ColorNode, ColorReplayError, ColorRowVerification, ColorRowVerificationError,
    ColorStructureRejection, ColorWriteError, MAX_COLOR_PARENTS, MAX_VALIDITY_AXES,
    MAX_WAIVER_CLOSURE_BYTES, MAX_WAIVER_DEPENDENCIES, NoSourceOriginVerifier, NoWaiverVerifier,
    PolicyDecision, SourceOrigin, SourceOriginRejection, SourceOriginRequest, SourceOriginVerifier,
    WAIVER_SCOPE_COLOR_UPGRADE, WAIVER_SCOPE_SOURCE_COLOR, Waiver, WaiverDependency, WaiverGrant,
    WaiverRejection, WaiverVerifier, verify_color_row_stream,
};
pub use crosswalk::{MAX_QTY_CROSSWALK_JSON_BYTES, QtyDimensionCrosswalkWrite};
pub use hash::{Blake3, ContentHash, hash_bytes};
pub use identity_migration::{
    ARTIFACT_CONTENT_IDENTITY_ROW_VERSION, ArtifactContentIdentity, ArtifactSemanticBinding,
    ArtifactSemanticBindingCandidates, ArtifactSemanticBindingWrite,
    EDGE_CONTENT_IDENTITY_ROW_VERSION, EdgeContentIdentity, EvidenceSemanticBinding,
    EvidenceSemanticBindingCandidates, EvidenceSemanticBindingWrite,
    IDENTITY_MIGRATION_RECEIPT_DOMAIN, IDENTITY_MIGRATION_RECEIPT_SCHEMA_DECLARATION,
    IDENTITY_MIGRATION_RECEIPT_VERSION, IDENTITY_MIGRATION_RECEIPT_WIRE_VERSION,
    IdentityMigrationCandidates, IdentityMigrationClaim, IdentityMigrationReceipt,
    IdentityMigrationReceiptId, IdentityMigrationReceiptSchemaV1, IdentityMigrationWireError,
    IdentityMigrationWrite, MAX_ARTIFACT_SEMANTIC_BINDING_CANDIDATES,
    MAX_EVIDENCE_SEMANTIC_BINDING_CANDIDATES, MAX_EVIDENCE_SEMANTIC_BINDING_NAME_BYTES,
    MAX_IDENTITY_MIGRATION_CANDIDATES, MAX_IDENTITY_MIGRATION_CONTEXT_BYTES,
    MAX_IDENTITY_MIGRATION_DOMAIN_BYTES, MAX_IDENTITY_MIGRATION_PAYLOAD_BYTES,
    MAX_IDENTITY_MIGRATION_RECEIPT_WIRE_BYTES, MAX_IDENTITY_MIGRATION_RULE_BYTES,
    MAX_IDENTITY_MIGRATION_SCHEMA_NAME_BYTES, OP_CONTENT_IDENTITY_ROW_VERSION, OpContentIdentity,
    TUNE_CONTENT_IDENTITY_ROW_VERSION, TuneContentIdentity,
};
pub use schema::{ALL_TABLES, SCHEMA_VERSION, STORAGE_CHUNK_LEN, V1_TABLES};
pub use state_checkpoint::{
    KnownStateSemantics, MAX_RUNTIME_STATE_CHECKPOINT_BYTES, MAX_STATE_CHECKPOINT_LAW_ID_BYTES,
    RUNTIME_STATE_ARTIFACT_KIND, STATE_CHECKPOINT_RECEIPT_IDENTITY_DOMAIN,
    STATE_CHECKPOINT_RECEIPT_IDENTITY_SCHEMA_DECLARATION,
    STATE_CHECKPOINT_RECEIPT_IDENTITY_VERSION, StateCheckpointClaim, StateCheckpointDecodeError,
    StateCheckpointReceipt, StateSlotId, VerifiedStateCheckpoint,
};
pub use travel::{
    BranchDiff, BranchInfo, ExecMode, ExplainNode, ExplainOp, GcReport, MAIN_BRANCH,
    ReplayMismatch, ReplayVerdict, ViewSnapshot,
};

use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

use fsqlite::{Connection, FrankenError, SqliteValue};

/// Maximum UTF-8 byte length of an artifact kind admitted or materialized.
pub const MAX_ARTIFACT_KIND_BYTES: usize = 256;

/// Maximum UTF-8 byte length of optional artifact metadata JSON.
pub const MAX_ARTIFACT_META_BYTES: usize = 1024 * 1024;

/// Conservative maximum byte length shared by each variable-size operation
/// field. Field-specific aliases keep the public contract independently
/// tunable without changing callers.
pub const MAX_OP_FIELD_BYTES: usize = 1024 * 1024;

/// Maximum byte length of an optional operation session identity.
pub const MAX_OP_SESSION_BYTES: usize = MAX_OP_FIELD_BYTES;

/// Maximum UTF-8 byte length of one operation IR JSON document.
pub const MAX_OP_IR_BYTES: usize = MAX_OP_FIELD_BYTES;

/// Maximum byte length of one operation RNG seed.
pub const MAX_OP_SEED_BYTES: usize = MAX_OP_FIELD_BYTES;

/// Maximum UTF-8 byte length of one operation versions JSON document.
pub const MAX_OP_VERSIONS_BYTES: usize = MAX_OP_FIELD_BYTES;

/// Maximum UTF-8 byte length of one operation budget JSON document.
pub const MAX_OP_BUDGET_BYTES: usize = MAX_OP_FIELD_BYTES;

/// Maximum UTF-8 byte length of one operation capability JSON document.
pub const MAX_OP_CAPABILITY_BYTES: usize = MAX_OP_FIELD_BYTES;

/// Maximum UTF-8 byte length of an optional operation diagnostic JSON document.
pub const MAX_OP_DIAG_BYTES: usize = MAX_OP_FIELD_BYTES;

/// Maximum byte length of a canonical autotuner kernel identity.
pub const MAX_TUNE_KERNEL_BYTES: usize = 64 * 1024;

/// Maximum byte length of a canonical autotuner shape-class identity.
pub const MAX_TUNE_SHAPE_CLASS_BYTES: usize = 64 * 1024;

/// Maximum byte length of an opaque autotuner machine fingerprint.
pub const MAX_TUNE_MACHINE_BYTES: usize = 256;

/// Maximum UTF-8 byte length of autotuner parameter JSON.
pub const MAX_TUNE_PARAMS_BYTES: usize = 1024 * 1024;

/// Maximum UTF-8 byte length of autotuner measurement JSON.
pub const MAX_TUNE_MEASURED_BYTES: usize = 1024 * 1024;

/// Maximum rows returned by one [`Ledger::tune_rows`] scan.
pub const MAX_TUNE_ROWS_PER_KERNEL: usize = 1024;

/// Maximum aggregate bytes returned by one [`Ledger::tune_rows`] scan.
pub const MAX_TUNE_SCAN_BYTES: usize = 16 * 1024 * 1024;

/// Maximum caller cap accepted by one bounded lineage-row query.
///
/// Each bounded query reads at most `cap + 1` fixed-size rows so it can report
/// truncation without materializing an unbounded producer or edge fan-out.
pub const MAX_LINEAGE_QUERY_ROWS: usize = 1024;

/// Crate version (compile-time stamp).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

// ---------------------------------------------------------------------------
// Semantic identity ownership (bead frankensim-semantic-identity-coverage-iu5l)
// ---------------------------------------------------------------------------

/// Schema version of the opaque physical-ledger identity.
pub const PHYSICAL_INSTANCE_IDENTITY_VERSION: u32 = 1;
/// Registry domain of the opaque physical-ledger identity.
pub const PHYSICAL_INSTANCE_IDENTITY_DOMAIN: &str = "org.frankensim.fs-ledger.physical-instance.v1";
/// Schema version of plain content-addressed artifact identity.
pub const ARTIFACT_CONTENT_IDENTITY_VERSION: u32 = 1;
/// Registry domain of plain content-addressed artifact identity.
pub const ARTIFACT_CONTENT_IDENTITY_DOMAIN: &str = "org.frankensim.fs-ledger.artifact-content.v1";
/// Schema version of the immutable session-mutation claim hash.
pub const SESSION_MUTATION_CLAIM_IDENTITY_VERSION: u32 = 1;
/// Registry domain of the immutable session-mutation claim hash.
pub const SESSION_MUTATION_CLAIM_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-ledger.session-mutation-claim.v1";
/// Schema version of the ordered session-terminal event hash.
pub const SESSION_TERMINAL_EVENTS_IDENTITY_VERSION: u32 = 2;
/// Registry domain of the ordered session-terminal event hash.
pub const SESSION_TERMINAL_EVENTS_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-ledger.session-terminal-events.v2";
/// Schema version of the complete session flush-batch witness.
pub const SESSION_FLUSH_BATCH_IDENTITY_VERSION: u32 = 2;
/// Registry domain of the complete session flush-batch witness.
pub const SESSION_FLUSH_BATCH_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-ledger.session-flush-batch.v2";
/// Schema version of a source-origin admission request.
pub const SOURCE_ORIGIN_REQUEST_IDENTITY_VERSION: u32 = 1;
/// Logical registry domain of a source-origin admission request.
pub const SOURCE_ORIGIN_REQUEST_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-ledger.source-origin-request.v1";
/// Schema version of a derived-color waiver signing subject.
pub const DERIVED_COLOR_WAIVER_SUBJECT_IDENTITY_VERSION: u32 = 3;
/// Logical registry domain of a derived-color waiver signing subject.
pub const DERIVED_COLOR_WAIVER_SUBJECT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-ledger.derived-color-waiver-subject.v3";
/// Schema version of a source-color waiver signing subject.
pub const SOURCE_COLOR_WAIVER_SUBJECT_IDENTITY_VERSION: u32 = 4;
/// Logical registry domain of a source-color waiver signing subject.
pub const SOURCE_COLOR_WAIVER_SUBJECT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-ledger.source-color-waiver-subject.v4";
/// Schema version of the complete color-node provenance hash.
pub const COLOR_NODE_IDENTITY_VERSION: u32 = 9;
/// Logical registry domain of the complete color-node provenance hash.
pub const COLOR_NODE_IDENTITY_DOMAIN: &str = "org.frankensim.fs-ledger.color-node.v9";
/// Schema version of the color-admission policy fingerprint.
pub const COLOR_ADMISSION_POLICY_IDENTITY_VERSION: u32 = 1;
/// Registry domain of the color-admission policy fingerprint.
pub const COLOR_ADMISSION_POLICY_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-ledger.color-admission-policy.v1";
/// Schema version of the persisted VCS ledger-lineage identity.
pub const VCS_LEDGER_LINEAGE_IDENTITY_VERSION: u32 = 1;
/// Registry domain of the persisted VCS ledger-lineage identity.
pub const VCS_LEDGER_LINEAGE_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-ledger.vcs-ledger-lineage.v1";
/// Schema version of a VCS semantic commit leaf.
pub const VCS_COMMIT_LEAF_IDENTITY_VERSION: u32 = 2;
/// Registry domain of a VCS semantic commit leaf.
pub const VCS_COMMIT_LEAF_IDENTITY_DOMAIN: &str = "org.frankensim.fs-ledger.vcs-commit-leaf.v2";
/// Schema version of a VCS semantic commit root.
pub const VCS_COMMIT_ROOT_IDENTITY_VERSION: u32 = 2;
/// Registry domain of a VCS semantic commit root.
pub const VCS_COMMIT_ROOT_IDENTITY_DOMAIN: &str = "org.frankensim.fs-ledger.vcs-commit-root.v2";
/// Schema version of the ledger/branch/root commit envelope key.
pub const VCS_COMMIT_ENVELOPE_IDENTITY_VERSION: u32 = 1;
/// Registry domain of the ledger/branch/root commit envelope key.
pub const VCS_COMMIT_ENVELOPE_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-ledger.vcs-commit-envelope.v1";

const SOURCE_ORIGIN_REQUEST_PREIMAGE_DOMAIN: &[u8] = b"frankensim/fs-ledger/source-origin-request";
const COLOR_WAIVER_SUBJECT_PREIMAGE_DOMAIN: &[u8] = b"frankensim/fs-ledger/color-waiver";
const COLOR_NODE_PREIMAGE_DOMAIN: &[u8] = b"frankensim/fs-ledger/color-node/v2";
const COLOR_ADMISSION_POLICY_PREIMAGE_DOMAIN: &str = "fs-ledger/color-admission-policy/v1";
const VCS_LEDGER_LINEAGE_PREIMAGE_DOMAIN: &[u8] = b"frankensim.fs-ledger.vcs.ledger-identity.v1";
const VCS_COMMIT_LEAF_PREIMAGE_DOMAIN: &[u8] = b"frankensim.fs-ledger.vcs.commit-leaf.v2";
const VCS_MERKLE_PAIR_PREIMAGE_DOMAIN: &[u8] = b"frankensim.fs-ledger.vcs.merkle-pair.v2";
const VCS_MERKLE_ODD_PREIMAGE_DOMAIN: &[u8] = b"frankensim.fs-ledger.vcs.merkle-odd.v2";
const VCS_COMMIT_ROOT_PREIMAGE_DOMAIN: &[u8] = b"frankensim.fs-ledger.vcs.commit-root.v2";

// These private witness shapes make every encoded input explicit in the owner
// file. Production encoders remain the authorities; declarations fingerprint
// those functions below, and the witnesses let the policy gate reject an
// unclassified field before generated registry data can drift.
#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct PhysicalInstanceIdentitySource {
    uuid: [u8; 16],
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct ArtifactContentIdentitySource {
    content: Vec<u8>,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct SessionMutationClaimIdentitySource {
    authority: ContentHash,
    ledger_instance_id: [u8; 16],
    governor_hash: ContentHash,
    session_open_hash: ContentHash,
    registry_schema_version: i64,
    kind: Vec<u8>,
    session: u64,
    ledger_scope: Vec<u8>,
    generation: u64,
    causal_ordinal: Option<u64>,
    payload: Vec<u8>,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct SessionTerminalEventIdentitySource {
    session: Vec<u8>,
    timestamp: i64,
    kind: Vec<u8>,
    payload: Option<Vec<u8>>,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct SessionFlushTerminalIdentitySource {
    authority: ContentHash,
    claim_hash: ContentHash,
    receipt_hash: ContentHash,
    event_count: usize,
    events_hash: ContentHash,
    encoded_bytes: usize,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct SessionFlushBatchIdentitySource {
    ledger_instance_id: [u8; 16],
    registry_schema_version: i64,
    terminals: Vec<SessionFlushTerminalIdentitySource>,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct SourceOriginRequestIdentitySource {
    node_name: Vec<u8>,
    claimed_color: Vec<u8>,
    origin: Vec<u8>,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct DerivedColorWaiverSubjectIdentitySource {
    operation_tag: u8,
    key_id: Vec<u8>,
    scope: Vec<u8>,
    node_name: Vec<u8>,
    claimed_color: Vec<u8>,
    annotation_id: Vec<u8>,
    annotation_signer: Vec<u8>,
    annotation_reason: Vec<u8>,
    parent_hashes: Vec<ContentHash>,
    expires_day: u32,
    signature: Vec<u8>,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct SourceColorWaiverSubjectIdentitySource {
    key_id: Vec<u8>,
    scope: Vec<u8>,
    node_name: Vec<u8>,
    claimed_color: Vec<u8>,
    annotation_id: Vec<u8>,
    annotation_signer: Vec<u8>,
    annotation_reason: Vec<u8>,
    parent_hashes: Vec<ContentHash>,
    expires_day: u32,
    signature: Vec<u8>,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct ColorNodeIdentitySource {
    node_id: u64,
    operation_tag: Option<u8>,
    name: Vec<u8>,
    color: Vec<u8>,
    parent_local_ids: Vec<u64>,
    parent_hashes: Vec<ContentHash>,
    demotions: Vec<Vec<u8>>,
    origin: Option<Vec<u8>>,
    origin_policy_fingerprint: Option<ContentHash>,
    waiver_dependencies: Vec<Vec<u8>>,
    waiver: Option<Vec<u8>>,
    grant_payload: Option<Vec<u8>>,
    grant_signature: Option<Vec<u8>>,
    waiver_policy_fingerprint: Option<ContentHash>,
    waiver_admission_day: Option<u32>,
    stored_hash: ContentHash,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct ColorAdmissionPolicyIdentitySource {
    color_write_row_schema_version: u32,
    color_algebra_version: u32,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct VcsLedgerLineageIdentitySource {
    mint_path: Vec<u8>,
    minted_ns: i64,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct VcsCommitEdgeIdentitySource {
    role: Vec<u8>,
    artifact_hash: ContentHash,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct VcsCommitLeafIdentitySource {
    ir: Vec<u8>,
    seed: Vec<u8>,
    versions: Vec<u8>,
    budget: Vec<u8>,
    capability: Vec<u8>,
    outcome: Option<Vec<u8>>,
    diagnostic: Option<Vec<u8>>,
    execution_mode: Vec<u8>,
    edges: Vec<VcsCommitEdgeIdentitySource>,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct VcsCommitRootIdentitySource {
    leaves: Vec<ContentHash>,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct VcsCommitEnvelopeIdentitySource {
    ledger: ContentHash,
    branch: i64,
    root: ContentHash,
}

#[allow(dead_code)]
fn identity_push_len(out: &mut Vec<u8>, len: usize) {
    out.extend_from_slice(
        &u64::try_from(len)
            .expect("identity witness length fits u64")
            .to_le_bytes(),
    );
}

#[allow(dead_code)]
fn identity_push_field(out: &mut Vec<u8>, bytes: &[u8]) {
    identity_push_len(out, bytes.len());
    out.extend_from_slice(bytes);
}

#[allow(dead_code)]
fn identity_update_len(hasher: &mut Blake3, len: usize) {
    hasher.update(
        &u64::try_from(len)
            .expect("identity witness length fits u64")
            .to_le_bytes(),
    );
}

#[allow(dead_code)]
fn identity_update_bytes(hasher: &mut Blake3, bytes: &[u8]) {
    identity_update_len(hasher, bytes.len());
    hasher.update(bytes);
}

#[allow(dead_code)]
fn identity_update_optional_u64(hasher: &mut Blake3, value: Option<u64>) {
    match value {
        Some(value) => {
            hasher.update(&[1]);
            hasher.update(&value.to_be_bytes());
        }
        None => hasher.update(&[0]),
    }
}

#[allow(dead_code)]
fn identity_vcs_hash_frame(hasher: &mut Blake3, bytes: &[u8]) {
    identity_update_len(hasher, bytes.len());
    hasher.update(bytes);
}

#[allow(dead_code)]
fn identity_vcs_domain_hasher(domain: &[u8]) -> Blake3 {
    let mut hasher = Blake3::new();
    identity_vcs_hash_frame(&mut hasher, b"domain");
    identity_vcs_hash_frame(&mut hasher, domain);
    hasher
}

#[allow(dead_code)]
fn identity_vcs_hash_field(hasher: &mut Blake3, name: &[u8], value: &[u8]) {
    identity_vcs_hash_frame(hasher, name);
    identity_vcs_hash_frame(hasher, value);
}

#[allow(dead_code)]
fn identity_vcs_hash_optional_field(hasher: &mut Blake3, name: &[u8], value: Option<&[u8]>) {
    identity_vcs_hash_frame(hasher, name);
    match value {
        Some(value) => {
            identity_vcs_hash_frame(hasher, b"present");
            identity_vcs_hash_frame(hasher, value);
        }
        None => identity_vcs_hash_frame(hasher, b"absent"),
    }
}

#[allow(dead_code)]
fn identity_vcs_framed_hash(domain: &[u8], fields: &[(&[u8], &[u8])]) -> ContentHash {
    let mut hasher = identity_vcs_domain_hasher(domain);
    for (name, value) in fields {
        identity_vcs_hash_field(&mut hasher, name, value);
    }
    hasher.finalize()
}

#[allow(dead_code)]
fn identity_schema_is_current(
    found_version: u32,
    found_domain: &str,
    supported_version: u32,
    supported_domain: &str,
) -> bool {
    found_version == supported_version && found_domain == supported_domain
}

#[allow(dead_code)]
fn ledger_physical_instance_identity(source: &PhysicalInstanceIdentitySource) -> [u8; 16] {
    source.uuid
}

#[allow(dead_code)]
fn ledger_artifact_content_identity(source: &ArtifactContentIdentitySource) -> ContentHash {
    hash_bytes(&source.content)
}

#[allow(dead_code)]
fn ledger_session_mutation_claim_identity_with_domain(
    source: &SessionMutationClaimIdentitySource,
    domain: &[u8],
) -> ContentHash {
    let mut hasher = Blake3::new();
    hasher.update(domain);
    hasher.update(source.authority.as_bytes());
    hasher.update(&source.ledger_instance_id);
    hasher.update(source.governor_hash.as_bytes());
    hasher.update(source.session_open_hash.as_bytes());
    hasher.update(&source.registry_schema_version.to_le_bytes());
    identity_update_bytes(&mut hasher, &source.kind);
    hasher.update(&source.session.to_be_bytes());
    identity_update_bytes(&mut hasher, &source.ledger_scope);
    hasher.update(&source.generation.to_be_bytes());
    identity_update_optional_u64(&mut hasher, source.causal_ordinal);
    hasher.update(hash_bytes(&source.payload).as_bytes());
    hasher.finalize()
}

#[allow(dead_code)]
fn ledger_session_mutation_claim_identity(
    source: &SessionMutationClaimIdentitySource,
) -> ContentHash {
    ledger_session_mutation_claim_identity_with_domain(
        source,
        b"org.frankensim.fs-ledger.session-mutation-claim.v1\0",
    )
}

#[allow(dead_code)]
fn ledger_session_terminal_events_identity_with_schema(
    events: &[SessionTerminalEventIdentitySource],
    declared_count: usize,
    domain: &[u8],
) -> ContentHash {
    let mut hasher = Blake3::new();
    hasher.update(domain);
    identity_update_len(&mut hasher, declared_count);
    for event in events {
        identity_update_bytes(&mut hasher, &event.session);
        hasher.update(&event.timestamp.to_le_bytes());
        identity_update_bytes(&mut hasher, &event.kind);
        match &event.payload {
            Some(payload) => {
                hasher.update(&[1]);
                identity_update_bytes(&mut hasher, payload);
            }
            None => hasher.update(&[0]),
        }
    }
    hasher.finalize()
}

#[allow(dead_code)]
fn ledger_session_terminal_events_identity(
    events: &[SessionTerminalEventIdentitySource],
) -> ContentHash {
    ledger_session_terminal_events_identity_with_schema(
        events,
        events.len(),
        b"org.frankensim.fs-ledger.session-terminal-events.v2\0",
    )
}

#[allow(dead_code)]
fn ledger_session_flush_batch_identity_with_schema(
    source: &SessionFlushBatchIdentitySource,
    declared_count: usize,
    domain: &[u8],
) -> ContentHash {
    let mut hasher = Blake3::new();
    hasher.update(domain);
    hasher.update(&source.ledger_instance_id);
    hasher.update(&source.registry_schema_version.to_le_bytes());
    identity_update_len(&mut hasher, declared_count);
    for terminal in &source.terminals {
        hasher.update(terminal.authority.as_bytes());
        hasher.update(terminal.claim_hash.as_bytes());
        hasher.update(terminal.receipt_hash.as_bytes());
        identity_update_len(&mut hasher, terminal.event_count);
        hasher.update(terminal.events_hash.as_bytes());
        identity_update_len(&mut hasher, terminal.encoded_bytes);
    }
    hasher.finalize()
}

#[allow(dead_code)]
fn ledger_session_flush_batch_identity(source: &SessionFlushBatchIdentitySource) -> ContentHash {
    ledger_session_flush_batch_identity_with_schema(
        source,
        source.terminals.len(),
        b"org.frankensim.fs-ledger.session-flush-batch.v2\0",
    )
}

#[allow(dead_code)]
fn ledger_source_origin_request_identity_with_schema(
    source: &SourceOriginRequestIdentitySource,
    version: u8,
    domain: &[u8],
) -> Vec<u8> {
    let mut out = vec![version];
    identity_push_field(&mut out, domain);
    identity_push_field(&mut out, &source.node_name);
    identity_push_field(&mut out, &source.claimed_color);
    out.extend_from_slice(&source.origin);
    out
}

#[allow(dead_code)]
fn ledger_source_origin_request_identity(source: &SourceOriginRequestIdentitySource) -> Vec<u8> {
    ledger_source_origin_request_identity_with_schema(
        source,
        SOURCE_ORIGIN_REQUEST_IDENTITY_VERSION as u8,
        SOURCE_ORIGIN_REQUEST_PREIMAGE_DOMAIN,
    )
}

#[allow(dead_code)]
fn identity_transport_has_versioned_domain(bytes: &[u8], version: u8, domain: &[u8]) -> bool {
    let Some((&found_version, rest)) = bytes.split_first() else {
        return false;
    };
    let Some(length_bytes) = rest.get(..8) else {
        return false;
    };
    let Ok(length_bytes) = <[u8; 8]>::try_from(length_bytes) else {
        return false;
    };
    let Ok(length) = usize::try_from(u64::from_le_bytes(length_bytes)) else {
        return false;
    };
    let Some(domain_end) = 8_usize.checked_add(length) else {
        return false;
    };
    found_version == version && length == domain.len() && rest.get(8..domain_end) == Some(domain)
}

#[allow(dead_code)]
fn ledger_source_origin_request_transport_guard(bytes: &[u8]) -> bool {
    identity_transport_has_versioned_domain(
        bytes,
        SOURCE_ORIGIN_REQUEST_IDENTITY_VERSION as u8,
        SOURCE_ORIGIN_REQUEST_PREIMAGE_DOMAIN,
    )
}

#[allow(dead_code)]
fn ledger_derived_color_waiver_subject_identity_with_schema(
    source: &DerivedColorWaiverSubjectIdentitySource,
    version: u8,
    domain: &[u8],
) -> Vec<u8> {
    let mut out = vec![version];
    identity_push_field(&mut out, domain);
    out.push(source.operation_tag);
    for field in [&source.key_id, &source.scope, &source.node_name] {
        identity_push_field(&mut out, field);
    }
    identity_push_field(&mut out, &source.claimed_color);
    for field in [
        &source.annotation_id,
        &source.annotation_signer,
        &source.annotation_reason,
    ] {
        identity_push_field(&mut out, field);
    }
    identity_push_len(&mut out, source.parent_hashes.len());
    for parent in &source.parent_hashes {
        out.extend_from_slice(parent.as_bytes());
    }
    out.extend_from_slice(&source.expires_day.to_le_bytes());
    out
}

#[allow(dead_code)]
fn ledger_derived_color_waiver_subject_identity(
    source: &DerivedColorWaiverSubjectIdentitySource,
) -> Vec<u8> {
    ledger_derived_color_waiver_subject_identity_with_schema(
        source,
        DERIVED_COLOR_WAIVER_SUBJECT_IDENTITY_VERSION as u8,
        COLOR_WAIVER_SUBJECT_PREIMAGE_DOMAIN,
    )
}

#[allow(dead_code)]
fn ledger_derived_color_waiver_subject_transport_guard(bytes: &[u8]) -> bool {
    identity_transport_has_versioned_domain(
        bytes,
        DERIVED_COLOR_WAIVER_SUBJECT_IDENTITY_VERSION as u8,
        COLOR_WAIVER_SUBJECT_PREIMAGE_DOMAIN,
    )
}

#[allow(dead_code)]
fn ledger_source_color_waiver_subject_identity_with_schema(
    source: &SourceColorWaiverSubjectIdentitySource,
    version: u8,
    domain: &[u8],
) -> Vec<u8> {
    let mut out = vec![version];
    identity_push_field(&mut out, domain);
    out.push(0);
    for field in [&source.key_id, &source.scope, &source.node_name] {
        identity_push_field(&mut out, field);
    }
    identity_push_field(&mut out, &source.claimed_color);
    for field in [
        &source.annotation_id,
        &source.annotation_signer,
        &source.annotation_reason,
    ] {
        identity_push_field(&mut out, field);
    }
    identity_push_len(&mut out, source.parent_hashes.len());
    for parent in &source.parent_hashes {
        out.extend_from_slice(parent.as_bytes());
    }
    out.extend_from_slice(&source.expires_day.to_le_bytes());
    out
}

#[allow(dead_code)]
fn ledger_source_color_waiver_subject_identity(
    source: &SourceColorWaiverSubjectIdentitySource,
) -> Vec<u8> {
    ledger_source_color_waiver_subject_identity_with_schema(
        source,
        SOURCE_COLOR_WAIVER_SUBJECT_IDENTITY_VERSION as u8,
        COLOR_WAIVER_SUBJECT_PREIMAGE_DOMAIN,
    )
}

#[allow(dead_code)]
fn ledger_source_color_waiver_subject_transport_guard(bytes: &[u8]) -> bool {
    identity_transport_has_versioned_domain(
        bytes,
        SOURCE_COLOR_WAIVER_SUBJECT_IDENTITY_VERSION as u8,
        COLOR_WAIVER_SUBJECT_PREIMAGE_DOMAIN,
    )
}

#[allow(dead_code)]
fn ledger_color_node_identity_with_schema(
    source: &ColorNodeIdentitySource,
    version: u8,
    domain: &[u8],
) -> ContentHash {
    let mut out = vec![version];
    identity_push_field(&mut out, domain);
    match source.operation_tag {
        Some(operation) => {
            out.push(1);
            out.push(operation);
        }
        None => out.push(0),
    }
    identity_push_field(&mut out, &source.name);
    identity_push_field(&mut out, &source.color);
    identity_push_len(&mut out, source.parent_hashes.len());
    for parent in &source.parent_hashes {
        identity_push_field(&mut out, parent.as_bytes());
    }
    identity_push_len(&mut out, source.demotions.len());
    for demotion in &source.demotions {
        out.extend_from_slice(demotion);
    }
    match &source.origin {
        Some(origin) => {
            out.push(1);
            out.extend_from_slice(origin);
        }
        None => out.push(0),
    }
    match source.origin_policy_fingerprint {
        Some(policy) => {
            out.push(1);
            out.extend_from_slice(policy.as_bytes());
        }
        None => out.push(0),
    }
    identity_push_len(&mut out, source.waiver_dependencies.len());
    for dependency in &source.waiver_dependencies {
        out.extend_from_slice(dependency);
    }
    match &source.waiver {
        Some(waiver) => {
            out.push(1);
            out.extend_from_slice(waiver);
        }
        None => out.push(0),
    }
    match (&source.grant_payload, &source.grant_signature) {
        (Some(payload), Some(signature)) => {
            out.push(1);
            identity_push_field(&mut out, payload);
            identity_push_field(&mut out, signature);
        }
        _ => out.push(0),
    }
    match source.waiver_policy_fingerprint {
        Some(policy) => {
            out.push(1);
            out.extend_from_slice(policy.as_bytes());
        }
        None => out.push(0),
    }
    match source.waiver_admission_day {
        Some(day) => {
            out.push(1);
            out.extend_from_slice(&day.to_le_bytes());
        }
        None => out.push(0),
    }
    hash_bytes(&out)
}

#[allow(dead_code)]
fn ledger_color_node_identity(source: &ColorNodeIdentitySource) -> ContentHash {
    ledger_color_node_identity_with_schema(
        source,
        COLOR_NODE_IDENTITY_VERSION as u8,
        COLOR_NODE_PREIMAGE_DOMAIN,
    )
}

#[allow(dead_code)]
fn ledger_color_admission_policy_identity_with_schema(
    source: &ColorAdmissionPolicyIdentitySource,
    domain: &str,
) -> ContentHash {
    hash_bytes(
        format!(
            "{domain}/row-schema={}/algebra={}",
            source.color_write_row_schema_version, source.color_algebra_version
        )
        .as_bytes(),
    )
}

#[allow(dead_code)]
fn ledger_color_admission_policy_identity(
    source: &ColorAdmissionPolicyIdentitySource,
) -> ContentHash {
    ledger_color_admission_policy_identity_with_schema(
        source,
        COLOR_ADMISSION_POLICY_PREIMAGE_DOMAIN,
    )
}

#[allow(dead_code)]
fn ledger_vcs_ledger_lineage_identity_with_domain(
    source: &VcsLedgerLineageIdentitySource,
    domain: &[u8],
) -> ContentHash {
    identity_vcs_framed_hash(
        domain,
        &[
            (b"path", &source.mint_path),
            (b"minted_ns", &source.minted_ns.to_le_bytes()),
        ],
    )
}

#[allow(dead_code)]
fn ledger_vcs_ledger_lineage_identity(source: &VcsLedgerLineageIdentitySource) -> ContentHash {
    ledger_vcs_ledger_lineage_identity_with_domain(source, VCS_LEDGER_LINEAGE_PREIMAGE_DOMAIN)
}

#[allow(dead_code)]
fn ledger_vcs_commit_leaf_identity_with_domain(
    source: &VcsCommitLeafIdentitySource,
    domain: &[u8],
) -> ContentHash {
    let mut hasher = identity_vcs_domain_hasher(domain);
    identity_vcs_hash_field(&mut hasher, b"ir", &source.ir);
    identity_vcs_hash_field(&mut hasher, b"seed", &source.seed);
    identity_vcs_hash_field(&mut hasher, b"versions", &source.versions);
    identity_vcs_hash_field(&mut hasher, b"budget", &source.budget);
    identity_vcs_hash_field(&mut hasher, b"capability", &source.capability);
    identity_vcs_hash_optional_field(&mut hasher, b"outcome", source.outcome.as_deref());
    identity_vcs_hash_optional_field(&mut hasher, b"diag", source.diagnostic.as_deref());
    identity_vcs_hash_field(&mut hasher, b"exec_mode", &source.execution_mode);
    identity_vcs_hash_field(
        &mut hasher,
        b"edge_count",
        &u64::try_from(source.edges.len())
            .expect("bounded edge count fits u64")
            .to_le_bytes(),
    );
    for edge in &source.edges {
        identity_vcs_hash_field(&mut hasher, b"edge_role", &edge.role);
        identity_vcs_hash_field(&mut hasher, b"artifact_hash", edge.artifact_hash.as_bytes());
    }
    hasher.finalize()
}

#[allow(dead_code)]
fn ledger_vcs_commit_leaf_identity(source: &VcsCommitLeafIdentitySource) -> ContentHash {
    ledger_vcs_commit_leaf_identity_with_domain(source, VCS_COMMIT_LEAF_PREIMAGE_DOMAIN)
}

#[allow(dead_code)]
fn ledger_vcs_commit_root_identity_with_domains(
    source: &VcsCommitRootIdentitySource,
    pair_domain: &[u8],
    odd_domain: &[u8],
    root_domain: &[u8],
) -> ContentHash {
    let leaf_count = u64::try_from(source.leaves.len())
        .expect("bounded leaf count fits u64")
        .to_le_bytes();
    if source.leaves.is_empty() {
        return identity_vcs_framed_hash(root_domain, &[(b"leaf_count", &leaf_count)]);
    }
    let mut level = source.leaves.clone();
    while level.len() > 1 {
        let mut next = Vec::with_capacity(level.len().div_ceil(2));
        for pair in level.chunks(2) {
            next.push(if pair.len() == 2 {
                identity_vcs_framed_hash(
                    pair_domain,
                    &[
                        (b"left", pair[0].as_bytes()),
                        (b"right", pair[1].as_bytes()),
                    ],
                )
            } else {
                identity_vcs_framed_hash(odd_domain, &[(b"child", pair[0].as_bytes())])
            });
        }
        level = next;
    }
    identity_vcs_framed_hash(
        root_domain,
        &[(b"leaf_count", &leaf_count), (b"tree", level[0].as_bytes())],
    )
}

#[allow(dead_code)]
fn ledger_vcs_commit_root_identity(source: &VcsCommitRootIdentitySource) -> ContentHash {
    ledger_vcs_commit_root_identity_with_domains(
        source,
        VCS_MERKLE_PAIR_PREIMAGE_DOMAIN,
        VCS_MERKLE_ODD_PREIMAGE_DOMAIN,
        VCS_COMMIT_ROOT_PREIMAGE_DOMAIN,
    )
}

#[allow(dead_code)]
fn ledger_vcs_commit_envelope_identity(source: &VcsCommitEnvelopeIdentitySource) -> Vec<u8> {
    let mut out = Vec::with_capacity(72);
    out.extend_from_slice(source.ledger.as_bytes());
    out.extend_from_slice(&source.branch.to_le_bytes());
    out.extend_from_slice(source.root.as_bytes());
    out
}

#[allow(dead_code)]
fn classify_physical_instance_identity_fields(source: &PhysicalInstanceIdentitySource) {
    let PhysicalInstanceIdentitySource { uuid } = source;
    let _ = uuid;
}

#[allow(dead_code)]
fn classify_artifact_content_identity_fields(source: &ArtifactContentIdentitySource) {
    let ArtifactContentIdentitySource { content } = source;
    let _ = content;
}

#[allow(dead_code)]
fn classify_session_mutation_claim_identity_fields(source: &SessionMutationClaimIdentitySource) {
    let SessionMutationClaimIdentitySource {
        authority,
        ledger_instance_id,
        governor_hash,
        session_open_hash,
        registry_schema_version,
        kind,
        session,
        ledger_scope,
        generation,
        causal_ordinal,
        payload,
    } = source;
    let _ = (
        authority,
        ledger_instance_id,
        governor_hash,
        session_open_hash,
        registry_schema_version,
        kind,
        session,
        ledger_scope,
        generation,
        causal_ordinal,
        payload,
    );
}

#[allow(dead_code)]
fn classify_session_terminal_events_identity_fields(source: &SessionTerminalEventIdentitySource) {
    let SessionTerminalEventIdentitySource {
        session,
        timestamp,
        kind,
        payload,
    } = source;
    let _ = (session, timestamp, kind, payload);
}

#[allow(dead_code)]
fn classify_session_flush_batch_identity_fields(
    batch: &SessionFlushBatchIdentitySource,
    terminal: &SessionFlushTerminalIdentitySource,
) {
    let SessionFlushBatchIdentitySource {
        ledger_instance_id,
        registry_schema_version,
        terminals,
    } = batch;
    let SessionFlushTerminalIdentitySource {
        authority,
        claim_hash,
        receipt_hash,
        event_count,
        events_hash,
        encoded_bytes,
    } = terminal;
    let _ = (
        ledger_instance_id,
        registry_schema_version,
        terminals,
        authority,
        claim_hash,
        receipt_hash,
        event_count,
        events_hash,
        encoded_bytes,
    );
}

#[allow(dead_code)]
fn classify_source_origin_request_identity_fields(source: &SourceOriginRequestIdentitySource) {
    let SourceOriginRequestIdentitySource {
        node_name,
        claimed_color,
        origin,
    } = source;
    let _ = (node_name, claimed_color, origin);
}

#[allow(dead_code)]
fn classify_derived_color_waiver_subject_identity_fields(
    source: &DerivedColorWaiverSubjectIdentitySource,
) {
    let DerivedColorWaiverSubjectIdentitySource {
        operation_tag,
        key_id,
        scope,
        node_name,
        claimed_color,
        annotation_id,
        annotation_signer,
        annotation_reason,
        parent_hashes,
        expires_day,
        signature,
    } = source;
    let _ = (
        operation_tag,
        key_id,
        scope,
        node_name,
        claimed_color,
        annotation_id,
        annotation_signer,
        annotation_reason,
        parent_hashes,
        expires_day,
        signature,
    );
}

#[allow(dead_code)]
fn classify_source_color_waiver_subject_identity_fields(
    source: &SourceColorWaiverSubjectIdentitySource,
) {
    let SourceColorWaiverSubjectIdentitySource {
        key_id,
        scope,
        node_name,
        claimed_color,
        annotation_id,
        annotation_signer,
        annotation_reason,
        parent_hashes,
        expires_day,
        signature,
    } = source;
    let _ = (
        key_id,
        scope,
        node_name,
        claimed_color,
        annotation_id,
        annotation_signer,
        annotation_reason,
        parent_hashes,
        expires_day,
        signature,
    );
}

#[allow(dead_code)]
fn classify_color_node_identity_fields(source: &ColorNodeIdentitySource) {
    let ColorNodeIdentitySource {
        node_id,
        operation_tag,
        name,
        color,
        parent_local_ids,
        parent_hashes,
        demotions,
        origin,
        origin_policy_fingerprint,
        waiver_dependencies,
        waiver,
        grant_payload,
        grant_signature,
        waiver_policy_fingerprint,
        waiver_admission_day,
        stored_hash,
    } = source;
    let _ = (
        node_id,
        operation_tag,
        name,
        color,
        parent_local_ids,
        parent_hashes,
        demotions,
        origin,
        origin_policy_fingerprint,
        waiver_dependencies,
        waiver,
        grant_payload,
        grant_signature,
        waiver_policy_fingerprint,
        waiver_admission_day,
        stored_hash,
    );
}

#[allow(dead_code)]
fn classify_color_admission_policy_identity_fields(source: &ColorAdmissionPolicyIdentitySource) {
    let ColorAdmissionPolicyIdentitySource {
        color_write_row_schema_version,
        color_algebra_version,
    } = source;
    let _ = (color_write_row_schema_version, color_algebra_version);
}

#[allow(dead_code)]
fn classify_vcs_ledger_lineage_identity_fields(source: &VcsLedgerLineageIdentitySource) {
    let VcsLedgerLineageIdentitySource {
        mint_path,
        minted_ns,
    } = source;
    let _ = (mint_path, minted_ns);
}

#[allow(dead_code)]
fn classify_vcs_commit_leaf_identity_fields(
    leaf: &VcsCommitLeafIdentitySource,
    edge: &VcsCommitEdgeIdentitySource,
) {
    let VcsCommitLeafIdentitySource {
        ir,
        seed,
        versions,
        budget,
        capability,
        outcome,
        diagnostic,
        execution_mode,
        edges,
    } = leaf;
    let VcsCommitEdgeIdentitySource {
        role,
        artifact_hash,
    } = edge;
    let _ = (
        ir,
        seed,
        versions,
        budget,
        capability,
        outcome,
        diagnostic,
        execution_mode,
        edges,
        role,
        artifact_hash,
    );
}

#[allow(dead_code)]
fn classify_vcs_commit_root_identity_fields(source: &VcsCommitRootIdentitySource) {
    let VcsCommitRootIdentitySource { leaves } = source;
    let _ = leaves;
}

#[allow(dead_code)]
fn classify_vcs_commit_envelope_identity_fields(source: &VcsCommitEnvelopeIdentitySource) {
    let VcsCommitEnvelopeIdentitySource {
        ledger,
        branch,
        root,
    } = source;
    let _ = (ledger, branch, root);
}

/// Owner declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const PHYSICAL_INSTANCE_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-ledger:physical-instance",
    "version_const=PHYSICAL_INSTANCE_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-ledger.physical-instance.v1",
    "domain_const=PHYSICAL_INSTANCE_IDENTITY_DOMAIN",
    "encoder=ledger_physical_instance_identity",
    "encoder_helpers=identity_schema_is_current",
    "schema_constants=PHYSICAL_INSTANCE_IDENTITY_VERSION,PHYSICAL_INSTANCE_IDENTITY_DOMAIN,crates/fs-ledger/src/schema.rs#V4,crates/fs-ledger/src/schema.rs#V5",
    "schema_functions=fresh_ledger_instance_id,decode_ledger_instance_id,normalize_schema_sql,LedgerInstanceId::as_bytes,Ledger::instance_id,Ledger::checked_instance_id,Ledger::attest_identity_guards,Ledger::open,Ledger::migrate,Ledger::migrate_from_observed_version,Ledger::seed_instance_id_if_missing,Ledger::read_current_instance_id,identity_schema_is_current",
    "schema_dependencies=none",
    "digest=none-opaque-rfc4122-uuid",
    "encoding=fixed-width-key",
    "sources=PhysicalInstanceIdentitySource",
    "source_fields=PhysicalInstanceIdentitySource.uuid:semantic",
    "source_bindings=PhysicalInstanceIdentitySource.uuid>uuid-bytes",
    "external_semantic_fields=none",
    "semantic_fields=uuid-bytes",
    "excluded_fields=path:location-is-not-instance-authority,path-alias:aliases-share-one-persisted-uuid,handle-address:object-location-is-ephemeral,reopen-count:reopen-is-an-envelope-event",
    "consumers=Ledger::instance_id,Ledger::checked_instance_id,Ledger::open,fs-ledger:session-mutation-claim,fs-ledger:session-flush-batch",
    "mutations=uuid-bytes:crates/fs-ledger/src/lib.rs#physical_instance_identity_fields_move_independently",
    "nonsemantic_mutations=path:crates/fs-ledger/src/lib.rs#physical_instance_excluded_fields_do_not_move_identity,path-alias:crates/fs-ledger/src/lib.rs#physical_instance_excluded_fields_do_not_move_identity,handle-address:crates/fs-ledger/src/lib.rs#physical_instance_excluded_fields_do_not_move_identity,reopen-count:crates/fs-ledger/src/lib.rs#physical_instance_excluded_fields_do_not_move_identity",
    "field_guard=classify_physical_instance_identity_fields",
    "transport_guard=ledger_physical_instance_identity",
    "version_guard=crates/fs-ledger/src/lib.rs#ledger_semantic_identity_versions_fail_closed",
    "coupling_surface=fs-ledger:physical-instance",
];

/// Owner declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const ARTIFACT_CONTENT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-ledger:artifact-content",
    "version_const=ARTIFACT_CONTENT_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-ledger.artifact-content.v1",
    "domain_const=ARTIFACT_CONTENT_IDENTITY_DOMAIN",
    "encoder=ledger_artifact_content_identity",
    "encoder_helpers=none",
    "schema_constants=ARTIFACT_CONTENT_IDENTITY_VERSION,ARTIFACT_CONTENT_IDENTITY_DOMAIN,crates/fs-ledger/src/identity_migration.rs#ARTIFACT_CONTENT_IDENTITY_ROW_VERSION,crates/fs-ledger/src/identity_migration.rs#EDGE_CONTENT_IDENTITY_ROW_VERSION,crates/fs-ledger/src/identity_migration.rs#OP_CONTENT_IDENTITY_ROW_VERSION,crates/fs-ledger/src/identity_migration.rs#TUNE_CONTENT_IDENTITY_ROW_VERSION,crates/fs-ledger/src/schema.rs#V14,crates/fs-ledger/src/schema.rs#V15,crates/fs-ledger/src/schema.rs#V16,crates/fs-ledger/src/schema.rs#V18,crates/fs-ledger/src/schema.rs#V19,crates/fs-blake3/src/lib.rs#IV,crates/fs-blake3/src/lib.rs#MSG_PERMUTATION,crates/fs-blake3/src/lib.rs#BLOCK_LEN,crates/fs-blake3/src/lib.rs#CHUNK_LEN,crates/fs-blake3/src/lib.rs#CHUNK_START,crates/fs-blake3/src/lib.rs#CHUNK_END,crates/fs-blake3/src/lib.rs#PARENT,crates/fs-blake3/src/lib.rs#ROOT,crates/fs-blake3/src/lib.rs#MAX_DEPTH",
    "schema_functions=crates/fs-blake3/src/lib.rs#hash_bytes,crates/fs-blake3/src/lib.rs#Blake3::new,crates/fs-blake3/src/lib.rs#Blake3::update,crates/fs-blake3/src/lib.rs#Blake3::finalize,crates/fs-blake3/src/identity.rs#ContentId::of_bytes,Ledger::put_artifact,Ledger::artifact_writer,ArtifactWriter::finish,ArtifactWriter::finish_inner,Ledger::insert_inline_artifact,Ledger::read_artifact_chunks_with_info,Ledger::tune_put,Ledger::tune_put_if_absent,Ledger::tune_get_inner,crates/fs-ledger/src/identity_migration.rs#derive_op_content_identity,crates/fs-ledger/src/identity_migration.rs#derive_tune_content_identity,crates/fs-ledger/src/identity_migration.rs#Ledger::artifact_content_identity,crates/fs-ledger/src/identity_migration.rs#Ledger::verify_artifact_content_identity_backfill,crates/fs-ledger/src/identity_migration.rs#Ledger::edge_content_identity,crates/fs-ledger/src/identity_migration.rs#Ledger::verify_edge_content_identity_backfill,crates/fs-ledger/src/identity_migration.rs#Ledger::write_op_content_identity,crates/fs-ledger/src/identity_migration.rs#Ledger::insert_op_content_identity,crates/fs-ledger/src/identity_migration.rs#Ledger::op_content_identity,crates/fs-ledger/src/identity_migration.rs#Ledger::reconcile_op_content_identity,crates/fs-ledger/src/identity_migration.rs#Ledger::backfill_and_verify_op_content_identities,crates/fs-ledger/src/identity_migration.rs#Ledger::write_tune_content_identity,crates/fs-ledger/src/identity_migration.rs#Ledger::upsert_tune_content_identity,crates/fs-ledger/src/identity_migration.rs#Ledger::tune_content_identity,crates/fs-ledger/src/identity_migration.rs#Ledger::reconcile_tune_content_identity,crates/fs-ledger/src/identity_migration.rs#Ledger::verify_tune_content_identity,crates/fs-ledger/src/identity_migration.rs#Ledger::backfill_and_verify_tune_content_identities,crates/fs-ledger/src/identity_migration.rs#Ledger::bind_artifact_semantic_identity,crates/fs-ledger/src/identity_migration.rs#Ledger::artifact_semantic_binding,crates/fs-ledger/src/identity_migration.rs#Ledger::verify_artifact_semantic_bindings,identity_schema_is_current",
    "schema_dependencies=none",
    "digest=blake3-256-plain-hash",
    "encoding=typed-binary",
    "sources=ArtifactContentIdentitySource",
    "source_fields=ArtifactContentIdentitySource.content:semantic",
    "source_bindings=ArtifactContentIdentitySource.content>content-bytes",
    "external_semantic_fields=none",
    "semantic_fields=content-bytes",
    "excluded_fields=kind:typed-envelope-not-content,metadata:provenance-envelope-not-content,created-at:wall-clock-envelope,chunk-boundaries:storage-layout-only",
    "consumers=Ledger::put_artifact,ArtifactWriter::finish,Ledger::get_artifact,Ledger::read_artifact_chunks,Ledger::verify_artifact_integrity,Ledger::artifact_content_identity,Ledger::edge_content_identity,Ledger::op_content_identity,Ledger::reconcile_op_content_identity,Ledger::tune_content_identity,Ledger::reconcile_tune_content_identity,Ledger::bind_artifact_semantic_identity,Ledger::artifact_semantic_binding,fs-ledger:vcs-commit-leaf",
    "mutations=content-bytes:crates/fs-ledger/src/lib.rs#artifact_content_identity_fields_move_independently",
    "nonsemantic_mutations=kind:crates/fs-ledger/src/lib.rs#artifact_content_excluded_fields_do_not_move_identity,metadata:crates/fs-ledger/src/lib.rs#artifact_content_excluded_fields_do_not_move_identity,created-at:crates/fs-ledger/src/lib.rs#artifact_content_excluded_fields_do_not_move_identity,chunk-boundaries:crates/fs-ledger/src/lib.rs#artifact_content_excluded_fields_do_not_move_identity",
    "field_guard=classify_artifact_content_identity_fields",
    "transport_guard=ledger_artifact_content_identity",
    "version_guard=crates/fs-ledger/src/lib.rs#ledger_semantic_identity_versions_fail_closed",
    "coupling_surface=fs-ledger:artifact-content",
];

/// Owner declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const SESSION_MUTATION_CLAIM_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-ledger:session-mutation-claim",
    "version_const=SESSION_MUTATION_CLAIM_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-ledger.session-mutation-claim.v1",
    "domain_const=SESSION_MUTATION_CLAIM_IDENTITY_DOMAIN",
    "encoder=ledger_session_mutation_claim_identity",
    "encoder_helpers=ledger_session_mutation_claim_identity_with_domain,identity_update_len,identity_update_bytes,identity_update_optional_u64",
    "schema_constants=SESSION_MUTATION_CLAIM_IDENTITY_VERSION,SESSION_MUTATION_CLAIM_IDENTITY_DOMAIN,crates/fs-ledger/src/session_registry.rs#SESSION_REGISTRY_ROW_SCHEMA_VERSION,crates/fs-ledger/src/session_registry.rs#SESSION_CLAIM_HASH_DOMAIN,crates/fs-ledger/src/session_registry.rs#MAX_SESSION_TERMINAL_KIND_BYTES,crates/fs-ledger/src/session_registry.rs#MAX_SESSION_TERMINAL_SCOPE_BYTES,crates/fs-ledger/src/session_registry.rs#MAX_SESSION_CLAIM_PAYLOAD_BYTES",
    "schema_functions=crates/fs-ledger/src/session_registry.rs#compute_claim_hash,crates/fs-ledger/src/session_registry.rs#require_bounded_ascii,crates/fs-ledger/src/session_registry.rs#validate_claim,crates/fs-ledger/src/session_registry.rs#decode_stored_session_claim,crates/fs-ledger/src/session_registry.rs#StoredSessionMutationClaim::matches,crates/fs-ledger/src/session_registry.rs#Ledger::claim_session_mutation,identity_schema_is_current",
    "schema_dependencies=fs-ledger:artifact-content,fs-ledger:physical-instance",
    "digest=blake3-256-domain-separated",
    "encoding=typed-binary",
    "sources=SessionMutationClaimIdentitySource",
    "source_fields=SessionMutationClaimIdentitySource.authority:semantic,SessionMutationClaimIdentitySource.ledger_instance_id:semantic,SessionMutationClaimIdentitySource.governor_hash:semantic,SessionMutationClaimIdentitySource.session_open_hash:semantic,SessionMutationClaimIdentitySource.registry_schema_version:semantic,SessionMutationClaimIdentitySource.kind:semantic,SessionMutationClaimIdentitySource.session:semantic,SessionMutationClaimIdentitySource.ledger_scope:semantic,SessionMutationClaimIdentitySource.generation:semantic,SessionMutationClaimIdentitySource.causal_ordinal:semantic,SessionMutationClaimIdentitySource.payload:semantic",
    "source_bindings=SessionMutationClaimIdentitySource.authority>authority,SessionMutationClaimIdentitySource.ledger_instance_id>ledger-instance-id,SessionMutationClaimIdentitySource.governor_hash>governor-hash,SessionMutationClaimIdentitySource.session_open_hash>session-open-hash,SessionMutationClaimIdentitySource.registry_schema_version>registry-schema-version,SessionMutationClaimIdentitySource.kind>kind-byte-count+kind-bytes,SessionMutationClaimIdentitySource.session>session,SessionMutationClaimIdentitySource.ledger_scope>ledger-scope-byte-count+ledger-scope-bytes,SessionMutationClaimIdentitySource.generation>generation,SessionMutationClaimIdentitySource.causal_ordinal>causal-ordinal-presence+causal-ordinal-value,SessionMutationClaimIdentitySource.payload>payload-bytes-via-blake3",
    "external_semantic_fields=identity-domain",
    "semantic_fields=identity-domain,authority,ledger-instance-id,governor-hash,session-open-hash,registry-schema-version,kind-byte-count,kind-bytes,session,ledger-scope-byte-count,ledger-scope-bytes,generation,causal-ordinal-presence,causal-ordinal-value,payload-bytes-via-blake3",
    "excluded_fields=claim-rowid:database-envelope-only,created-at:wall-clock-envelope,terminalization-permit:execution-authority-not-claim-content",
    "consumers=Ledger::claim_session_mutation,Ledger::pending_session_mutation,Ledger::session_mutation_claim,Ledger::append_session_terminal_batch",
    "mutations=identity-domain:crates/fs-ledger/src/lib.rs#session_mutation_claim_identity_fields_move_independently,authority:crates/fs-ledger/src/lib.rs#session_mutation_claim_identity_fields_move_independently,ledger-instance-id:crates/fs-ledger/src/lib.rs#session_mutation_claim_identity_fields_move_independently,governor-hash:crates/fs-ledger/src/lib.rs#session_mutation_claim_identity_fields_move_independently,session-open-hash:crates/fs-ledger/src/lib.rs#session_mutation_claim_identity_fields_move_independently,registry-schema-version:crates/fs-ledger/src/lib.rs#session_mutation_claim_identity_fields_move_independently,kind-byte-count:crates/fs-ledger/src/lib.rs#session_mutation_claim_identity_fields_move_independently,kind-bytes:crates/fs-ledger/src/lib.rs#session_mutation_claim_identity_fields_move_independently,session:crates/fs-ledger/src/lib.rs#session_mutation_claim_identity_fields_move_independently,ledger-scope-byte-count:crates/fs-ledger/src/lib.rs#session_mutation_claim_identity_fields_move_independently,ledger-scope-bytes:crates/fs-ledger/src/lib.rs#session_mutation_claim_identity_fields_move_independently,generation:crates/fs-ledger/src/lib.rs#session_mutation_claim_identity_fields_move_independently,causal-ordinal-presence:crates/fs-ledger/src/lib.rs#session_mutation_claim_identity_fields_move_independently,causal-ordinal-value:crates/fs-ledger/src/lib.rs#session_mutation_claim_identity_fields_move_independently,payload-bytes-via-blake3:crates/fs-ledger/src/lib.rs#session_mutation_claim_identity_fields_move_independently",
    "nonsemantic_mutations=claim-rowid:crates/fs-ledger/src/lib.rs#session_mutation_claim_excluded_fields_do_not_move_identity,created-at:crates/fs-ledger/src/lib.rs#session_mutation_claim_excluded_fields_do_not_move_identity,terminalization-permit:crates/fs-ledger/src/lib.rs#session_mutation_claim_excluded_fields_do_not_move_identity",
    "field_guard=classify_session_mutation_claim_identity_fields",
    "transport_guard=ledger_session_mutation_claim_identity",
    "version_guard=crates/fs-ledger/src/lib.rs#ledger_semantic_identity_versions_fail_closed",
    "coupling_surface=fs-ledger:session-mutation-claim",
];

/// Owner declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const SESSION_TERMINAL_EVENTS_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-ledger:session-terminal-events",
    "version_const=SESSION_TERMINAL_EVENTS_IDENTITY_VERSION",
    "version=2",
    "domain=org.frankensim.fs-ledger.session-terminal-events.v2",
    "domain_const=SESSION_TERMINAL_EVENTS_IDENTITY_DOMAIN",
    "encoder=ledger_session_terminal_events_identity",
    "encoder_helpers=ledger_session_terminal_events_identity_with_schema",
    "schema_constants=SESSION_TERMINAL_EVENTS_IDENTITY_VERSION,SESSION_TERMINAL_EVENTS_IDENTITY_DOMAIN,crates/fs-ledger/src/session_registry.rs#SESSION_EVENTS_HASH_DOMAIN,crates/fs-ledger/src/session_registry.rs#MAX_SESSION_TERMINAL_EVENT_KIND_BYTES,crates/fs-ledger/src/session_registry.rs#MAX_SESSION_TERMINAL_EVENT_PAYLOAD_BYTES,crates/fs-ledger/src/session_registry.rs#MAX_SESSION_FLUSH_EVENTS",
    "schema_functions=crates/fs-ledger/src/session_registry.rs#events_hasher,crates/fs-ledger/src/session_registry.rs#update_event_preimage,crates/fs-ledger/src/session_registry.rs#hash_events,crates/fs-ledger/src/session_registry.rs#session_terminal_events_hash,crates/fs-ledger/src/session_registry.rs#validate_event,identity_schema_is_current",
    "schema_dependencies=none",
    "digest=blake3-256-domain-separated",
    "encoding=typed-binary",
    "sources=SessionTerminalEventIdentitySource",
    "source_fields=SessionTerminalEventIdentitySource.session:semantic,SessionTerminalEventIdentitySource.timestamp:semantic,SessionTerminalEventIdentitySource.kind:semantic,SessionTerminalEventIdentitySource.payload:semantic",
    "source_bindings=SessionTerminalEventIdentitySource.session>session-byte-count+session-bytes,SessionTerminalEventIdentitySource.timestamp>timestamp,SessionTerminalEventIdentitySource.kind>kind-byte-count+kind-bytes,SessionTerminalEventIdentitySource.payload>payload-presence+payload-byte-count+payload-bytes",
    "external_semantic_fields=identity-domain,event-count,event-order",
    "semantic_fields=identity-domain,event-count,event-order,session-byte-count,session-bytes,timestamp,kind-byte-count,kind-bytes,payload-presence,payload-byte-count,payload-bytes",
    "excluded_fields=global-event-rowid:storage-envelope-only,owner-link-rowid:storage-envelope-only,batch-rowid:storage-envelope-only",
    "consumers=session_terminal_events_hash,Ledger::append_session_terminal_batch,Ledger::session_mutation_terminal",
    "mutations=identity-domain:crates/fs-ledger/src/lib.rs#session_terminal_events_identity_fields_move_independently,event-count:crates/fs-ledger/src/lib.rs#session_terminal_events_identity_fields_move_independently,event-order:crates/fs-ledger/src/lib.rs#session_terminal_events_identity_fields_move_independently,session-byte-count:crates/fs-ledger/src/lib.rs#session_terminal_events_identity_fields_move_independently,session-bytes:crates/fs-ledger/src/lib.rs#session_terminal_events_identity_fields_move_independently,timestamp:crates/fs-ledger/src/lib.rs#session_terminal_events_identity_fields_move_independently,kind-byte-count:crates/fs-ledger/src/lib.rs#session_terminal_events_identity_fields_move_independently,kind-bytes:crates/fs-ledger/src/lib.rs#session_terminal_events_identity_fields_move_independently,payload-presence:crates/fs-ledger/src/lib.rs#session_terminal_events_identity_fields_move_independently,payload-byte-count:crates/fs-ledger/src/lib.rs#session_terminal_events_identity_fields_move_independently,payload-bytes:crates/fs-ledger/src/lib.rs#session_terminal_events_identity_fields_move_independently",
    "nonsemantic_mutations=global-event-rowid:crates/fs-ledger/src/lib.rs#session_terminal_events_excluded_fields_do_not_move_identity,owner-link-rowid:crates/fs-ledger/src/lib.rs#session_terminal_events_excluded_fields_do_not_move_identity,batch-rowid:crates/fs-ledger/src/lib.rs#session_terminal_events_excluded_fields_do_not_move_identity",
    "field_guard=classify_session_terminal_events_identity_fields",
    "transport_guard=ledger_session_terminal_events_identity",
    "version_guard=crates/fs-ledger/src/lib.rs#ledger_semantic_identity_versions_fail_closed",
    "coupling_surface=fs-ledger:session-terminal-events",
];

/// Owner declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const SESSION_FLUSH_BATCH_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-ledger:session-flush-batch",
    "version_const=SESSION_FLUSH_BATCH_IDENTITY_VERSION",
    "version=2",
    "domain=org.frankensim.fs-ledger.session-flush-batch.v2",
    "domain_const=SESSION_FLUSH_BATCH_IDENTITY_DOMAIN",
    "encoder=ledger_session_flush_batch_identity",
    "encoder_helpers=ledger_session_flush_batch_identity_with_schema",
    "schema_constants=SESSION_FLUSH_BATCH_IDENTITY_VERSION,SESSION_FLUSH_BATCH_IDENTITY_DOMAIN,crates/fs-ledger/src/session_registry.rs#SESSION_BATCH_HASH_DOMAIN,crates/fs-ledger/src/session_registry.rs#SESSION_REGISTRY_ROW_SCHEMA_VERSION,crates/fs-ledger/src/session_registry.rs#MAX_SESSION_FLUSH_TERMINALS,crates/fs-ledger/src/session_registry.rs#MAX_SESSION_FLUSH_EVENTS,crates/fs-ledger/src/session_registry.rs#MAX_SESSION_FLUSH_ENCODED_BYTES",
    "schema_functions=crates/fs-ledger/src/session_registry.rs#prepare_batch,crates/fs-ledger/src/session_registry.rs#update_prepared_terminal_preimage,crates/fs-ledger/src/session_registry.rs#validate_terminal,crates/fs-ledger/src/session_registry.rs#StoredSessionTerminal::matches,crates/fs-ledger/src/session_registry.rs#StoredSessionFlushBatch::matches,crates/fs-ledger/src/session_registry.rs#Ledger::append_session_terminal_batch,crates/fs-ledger/src/session_registry.rs#Ledger::verify_session_flush_batch_members,identity_schema_is_current",
    "schema_dependencies=fs-ledger:artifact-content,fs-ledger:physical-instance,fs-ledger:session-mutation-claim,fs-ledger:session-terminal-events",
    "digest=blake3-256-domain-separated",
    "encoding=typed-binary",
    "sources=SessionFlushBatchIdentitySource,SessionFlushTerminalIdentitySource",
    "source_fields=SessionFlushBatchIdentitySource.ledger_instance_id:semantic,SessionFlushBatchIdentitySource.registry_schema_version:semantic,SessionFlushBatchIdentitySource.terminals:semantic,SessionFlushTerminalIdentitySource.authority:semantic,SessionFlushTerminalIdentitySource.claim_hash:semantic,SessionFlushTerminalIdentitySource.receipt_hash:semantic,SessionFlushTerminalIdentitySource.event_count:semantic,SessionFlushTerminalIdentitySource.events_hash:semantic,SessionFlushTerminalIdentitySource.encoded_bytes:semantic",
    "source_bindings=SessionFlushBatchIdentitySource.ledger_instance_id>ledger-instance-id,SessionFlushBatchIdentitySource.registry_schema_version>registry-schema-version,SessionFlushBatchIdentitySource.terminals>terminal-count+terminal-order,SessionFlushTerminalIdentitySource.authority>authority,SessionFlushTerminalIdentitySource.claim_hash>claim-hash,SessionFlushTerminalIdentitySource.receipt_hash>receipt-hash,SessionFlushTerminalIdentitySource.event_count>event-count,SessionFlushTerminalIdentitySource.events_hash>events-hash,SessionFlushTerminalIdentitySource.encoded_bytes>encoded-byte-count",
    "external_semantic_fields=identity-domain",
    "semantic_fields=identity-domain,ledger-instance-id,registry-schema-version,terminal-count,terminal-order,authority,claim-hash,receipt-hash,event-count,events-hash,encoded-byte-count",
    "excluded_fields=batch-rowid:storage-envelope-only,created-at:wall-clock-envelope,terminal-rowid:storage-envelope-only,terminalization-permit:execution-authority-not-batch-content",
    "consumers=Ledger::append_session_terminal_batch,Ledger::session_mutation_terminal,Ledger::session_flush_batch",
    "mutations=identity-domain:crates/fs-ledger/src/lib.rs#session_flush_batch_identity_fields_move_independently,ledger-instance-id:crates/fs-ledger/src/lib.rs#session_flush_batch_identity_fields_move_independently,registry-schema-version:crates/fs-ledger/src/lib.rs#session_flush_batch_identity_fields_move_independently,terminal-count:crates/fs-ledger/src/lib.rs#session_flush_batch_identity_fields_move_independently,terminal-order:crates/fs-ledger/src/lib.rs#session_flush_batch_identity_fields_move_independently,authority:crates/fs-ledger/src/lib.rs#session_flush_batch_identity_fields_move_independently,claim-hash:crates/fs-ledger/src/lib.rs#session_flush_batch_identity_fields_move_independently,receipt-hash:crates/fs-ledger/src/lib.rs#session_flush_batch_identity_fields_move_independently,event-count:crates/fs-ledger/src/lib.rs#session_flush_batch_identity_fields_move_independently,events-hash:crates/fs-ledger/src/lib.rs#session_flush_batch_identity_fields_move_independently,encoded-byte-count:crates/fs-ledger/src/lib.rs#session_flush_batch_identity_fields_move_independently",
    "nonsemantic_mutations=batch-rowid:crates/fs-ledger/src/lib.rs#session_flush_batch_excluded_fields_do_not_move_identity,created-at:crates/fs-ledger/src/lib.rs#session_flush_batch_excluded_fields_do_not_move_identity,terminal-rowid:crates/fs-ledger/src/lib.rs#session_flush_batch_excluded_fields_do_not_move_identity,terminalization-permit:crates/fs-ledger/src/lib.rs#session_flush_batch_excluded_fields_do_not_move_identity",
    "field_guard=classify_session_flush_batch_identity_fields",
    "transport_guard=ledger_session_flush_batch_identity",
    "version_guard=crates/fs-ledger/src/lib.rs#ledger_semantic_identity_versions_fail_closed",
    "coupling_surface=fs-ledger:session-flush-batch",
];

/// Owner declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const SOURCE_ORIGIN_REQUEST_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-ledger:source-origin-request",
    "version_const=SOURCE_ORIGIN_REQUEST_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-ledger.source-origin-request.v1",
    "domain_const=SOURCE_ORIGIN_REQUEST_IDENTITY_DOMAIN",
    "encoder=ledger_source_origin_request_identity",
    "encoder_helpers=ledger_source_origin_request_identity_with_schema,identity_push_len,identity_push_field,identity_transport_has_versioned_domain",
    "schema_constants=SOURCE_ORIGIN_REQUEST_IDENTITY_VERSION,SOURCE_ORIGIN_REQUEST_IDENTITY_DOMAIN,SOURCE_ORIGIN_REQUEST_PREIMAGE_DOMAIN,crates/fs-ledger/src/colors.rs#SOURCE_ORIGIN_REQUEST_DOMAIN,crates/fs-evidence/src/color.rs#COLOR_ALGEBRA_VERSION",
    "schema_functions=crates/fs-ledger/src/colors.rs#SourceOriginRequest::canonical_bytes,crates/fs-ledger/src/colors.rs#push_source_origin,crates/fs-ledger/src/colors.rs#source_origin_canonical_bytes,crates/fs-ledger/src/colors.rs#push_len,crates/fs-ledger/src/colors.rs#push_field,crates/fs-ledger/src/colors.rs#numerical_kind_tag,crates/fs-evidence/src/color.rs#Color::canonical_bytes,identity_transport_has_versioned_domain,identity_schema_is_current",
    "schema_dependencies=fs-ledger:artifact-content",
    "digest=none-exact-canonical-signing-transport",
    "encoding=canonical-transport-exact-bits",
    "sources=SourceOriginRequestIdentitySource",
    "source_fields=SourceOriginRequestIdentitySource.node_name:semantic,SourceOriginRequestIdentitySource.claimed_color:semantic,SourceOriginRequestIdentitySource.origin:semantic",
    "source_bindings=SourceOriginRequestIdentitySource.node_name>node-name-byte-count+node-name,SourceOriginRequestIdentitySource.claimed_color>claimed-color-byte-count+claimed-color-canonical-bytes,SourceOriginRequestIdentitySource.origin>typed-origin-canonical-bytes",
    "external_semantic_fields=transport-version,domain-byte-count,preimage-domain",
    "semantic_fields=transport-version,domain-byte-count,preimage-domain,node-name-byte-count,node-name,claimed-color-byte-count,claimed-color-canonical-bytes,typed-origin-canonical-bytes",
    "excluded_fields=verifier-policy-fingerprint:callback-result-not-request,callback-order:execution-envelope-only",
    "consumers=SourceOriginRequest::canonical_bytes,SourceOriginVerifier::verify,ColorGraph::source_with_origin",
    "mutations=transport-version:crates/fs-ledger/src/lib.rs#source_origin_request_identity_fields_move_independently,domain-byte-count:crates/fs-ledger/src/lib.rs#source_origin_request_identity_fields_move_independently,preimage-domain:crates/fs-ledger/src/lib.rs#source_origin_request_identity_fields_move_independently,node-name-byte-count:crates/fs-ledger/src/lib.rs#source_origin_request_identity_fields_move_independently,node-name:crates/fs-ledger/src/lib.rs#source_origin_request_identity_fields_move_independently,claimed-color-byte-count:crates/fs-ledger/src/lib.rs#source_origin_request_identity_fields_move_independently,claimed-color-canonical-bytes:crates/fs-ledger/src/lib.rs#source_origin_request_identity_fields_move_independently,typed-origin-canonical-bytes:crates/fs-ledger/src/lib.rs#source_origin_request_identity_fields_move_independently",
    "nonsemantic_mutations=verifier-policy-fingerprint:crates/fs-ledger/src/lib.rs#source_origin_request_excluded_fields_do_not_move_identity,callback-order:crates/fs-ledger/src/lib.rs#source_origin_request_excluded_fields_do_not_move_identity",
    "field_guard=classify_source_origin_request_identity_fields",
    "transport_guard=ledger_source_origin_request_transport_guard",
    "version_guard=crates/fs-ledger/src/lib.rs#ledger_semantic_identity_versions_fail_closed",
    "coupling_surface=fs-ledger:source-origin-request",
];

/// Owner declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const DERIVED_COLOR_WAIVER_SUBJECT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-ledger:derived-color-waiver-subject",
    "version_const=DERIVED_COLOR_WAIVER_SUBJECT_IDENTITY_VERSION",
    "version=3",
    "domain=org.frankensim.fs-ledger.derived-color-waiver-subject.v3",
    "domain_const=DERIVED_COLOR_WAIVER_SUBJECT_IDENTITY_DOMAIN",
    "encoder=ledger_derived_color_waiver_subject_identity",
    "encoder_helpers=ledger_derived_color_waiver_subject_identity_with_schema",
    "schema_constants=DERIVED_COLOR_WAIVER_SUBJECT_IDENTITY_VERSION,DERIVED_COLOR_WAIVER_SUBJECT_IDENTITY_DOMAIN,COLOR_WAIVER_SUBJECT_PREIMAGE_DOMAIN,crates/fs-ledger/src/colors.rs#WAIVER_PAYLOAD_DOMAIN,crates/fs-ledger/src/colors.rs#WAIVER_SCOPE_COLOR_UPGRADE",
    "schema_functions=crates/fs-ledger/src/colors.rs#WaiverGrant::signing_payload,crates/fs-ledger/src/colors.rs#WaiverGrant::signing_payload_for,crates/fs-ledger/src/colors.rs#WaiverGrant::payload_version,crates/fs-ledger/src/colors.rs#interval_op_tag,crates/fs-ledger/src/colors.rs#push_len,crates/fs-ledger/src/colors.rs#push_field,crates/fs-ledger/src/colors.rs#validate_waiver_grant,identity_transport_has_versioned_domain,identity_schema_is_current",
    "schema_dependencies=none",
    "digest=none-exact-canonical-signing-transport",
    "encoding=canonical-transport-exact-bits",
    "sources=DerivedColorWaiverSubjectIdentitySource",
    "source_fields=DerivedColorWaiverSubjectIdentitySource.operation_tag:semantic,DerivedColorWaiverSubjectIdentitySource.key_id:semantic,DerivedColorWaiverSubjectIdentitySource.scope:semantic,DerivedColorWaiverSubjectIdentitySource.node_name:semantic,DerivedColorWaiverSubjectIdentitySource.claimed_color:semantic,DerivedColorWaiverSubjectIdentitySource.annotation_id:semantic,DerivedColorWaiverSubjectIdentitySource.annotation_signer:semantic,DerivedColorWaiverSubjectIdentitySource.annotation_reason:semantic,DerivedColorWaiverSubjectIdentitySource.parent_hashes:semantic,DerivedColorWaiverSubjectIdentitySource.expires_day:semantic,DerivedColorWaiverSubjectIdentitySource.signature:nonsemantic:self-signature-is-outside-its-subject",
    "source_bindings=DerivedColorWaiverSubjectIdentitySource.operation_tag>operation-tag,DerivedColorWaiverSubjectIdentitySource.key_id>key-id-byte-count+key-id,DerivedColorWaiverSubjectIdentitySource.scope>scope-byte-count+scope,DerivedColorWaiverSubjectIdentitySource.node_name>node-name-byte-count+node-name,DerivedColorWaiverSubjectIdentitySource.claimed_color>claimed-color-byte-count+claimed-color-canonical-bytes,DerivedColorWaiverSubjectIdentitySource.annotation_id>annotation-id-byte-count+annotation-id,DerivedColorWaiverSubjectIdentitySource.annotation_signer>annotation-signer-byte-count+annotation-signer,DerivedColorWaiverSubjectIdentitySource.annotation_reason>annotation-reason-byte-count+annotation-reason,DerivedColorWaiverSubjectIdentitySource.parent_hashes>parent-count+parent-order+parent-hashes,DerivedColorWaiverSubjectIdentitySource.expires_day>expires-day",
    "external_semantic_fields=transport-version,domain-byte-count,preimage-domain",
    "semantic_fields=transport-version,domain-byte-count,preimage-domain,operation-tag,key-id-byte-count,key-id,scope-byte-count,scope,node-name-byte-count,node-name,claimed-color-byte-count,claimed-color-canonical-bytes,annotation-id-byte-count,annotation-id,annotation-signer-byte-count,annotation-signer,annotation-reason-byte-count,annotation-reason,parent-count,parent-order,parent-hashes,expires-day",
    "excluded_fields=admission-day:verification-context-not-signed-subject,policy-fingerprint:verifier-result-not-signed-subject",
    "consumers=WaiverGrant::signing_payload,WaiverVerifier::verify,ColorGraph::derive_waived,ColorGraph::verify_replay",
    "mutations=transport-version:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_identity_fields_move_independently,domain-byte-count:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_identity_fields_move_independently,preimage-domain:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_identity_fields_move_independently,operation-tag:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_identity_fields_move_independently,key-id-byte-count:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_identity_fields_move_independently,key-id:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_identity_fields_move_independently,scope-byte-count:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_identity_fields_move_independently,scope:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_identity_fields_move_independently,node-name-byte-count:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_identity_fields_move_independently,node-name:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_identity_fields_move_independently,claimed-color-byte-count:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_identity_fields_move_independently,claimed-color-canonical-bytes:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_identity_fields_move_independently,annotation-id-byte-count:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_identity_fields_move_independently,annotation-id:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_identity_fields_move_independently,annotation-signer-byte-count:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_identity_fields_move_independently,annotation-signer:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_identity_fields_move_independently,annotation-reason-byte-count:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_identity_fields_move_independently,annotation-reason:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_identity_fields_move_independently,parent-count:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_identity_fields_move_independently,parent-order:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_identity_fields_move_independently,parent-hashes:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_identity_fields_move_independently,expires-day:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_identity_fields_move_independently",
    "nonsemantic_mutations=DerivedColorWaiverSubjectIdentitySource.signature:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_excluded_fields_do_not_move_identity,admission-day:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_excluded_fields_do_not_move_identity,policy-fingerprint:crates/fs-ledger/src/lib.rs#derived_color_waiver_subject_excluded_fields_do_not_move_identity",
    "field_guard=classify_derived_color_waiver_subject_identity_fields",
    "transport_guard=ledger_derived_color_waiver_subject_transport_guard",
    "version_guard=crates/fs-ledger/src/lib.rs#ledger_semantic_identity_versions_fail_closed",
    "coupling_surface=fs-ledger:derived-color-waiver-subject",
];

/// Owner declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const SOURCE_COLOR_WAIVER_SUBJECT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-ledger:source-color-waiver-subject",
    "version_const=SOURCE_COLOR_WAIVER_SUBJECT_IDENTITY_VERSION",
    "version=4",
    "domain=org.frankensim.fs-ledger.source-color-waiver-subject.v4",
    "domain_const=SOURCE_COLOR_WAIVER_SUBJECT_IDENTITY_DOMAIN",
    "encoder=ledger_source_color_waiver_subject_identity",
    "encoder_helpers=ledger_source_color_waiver_subject_identity_with_schema",
    "schema_constants=SOURCE_COLOR_WAIVER_SUBJECT_IDENTITY_VERSION,SOURCE_COLOR_WAIVER_SUBJECT_IDENTITY_DOMAIN,COLOR_WAIVER_SUBJECT_PREIMAGE_DOMAIN,crates/fs-ledger/src/colors.rs#WAIVER_PAYLOAD_DOMAIN,crates/fs-ledger/src/colors.rs#WAIVER_SCOPE_SOURCE_COLOR",
    "schema_functions=crates/fs-ledger/src/colors.rs#WaiverGrant::signing_payload_source,crates/fs-ledger/src/colors.rs#WaiverGrant::signing_payload_for,crates/fs-ledger/src/colors.rs#WaiverGrant::payload_version,crates/fs-ledger/src/colors.rs#push_len,crates/fs-ledger/src/colors.rs#push_field,crates/fs-ledger/src/colors.rs#validate_waiver_grant,identity_transport_has_versioned_domain,identity_schema_is_current",
    "schema_dependencies=none",
    "digest=none-exact-canonical-signing-transport",
    "encoding=canonical-transport-exact-bits",
    "sources=SourceColorWaiverSubjectIdentitySource",
    "source_fields=SourceColorWaiverSubjectIdentitySource.key_id:semantic,SourceColorWaiverSubjectIdentitySource.scope:semantic,SourceColorWaiverSubjectIdentitySource.node_name:semantic,SourceColorWaiverSubjectIdentitySource.claimed_color:semantic,SourceColorWaiverSubjectIdentitySource.annotation_id:semantic,SourceColorWaiverSubjectIdentitySource.annotation_signer:semantic,SourceColorWaiverSubjectIdentitySource.annotation_reason:semantic,SourceColorWaiverSubjectIdentitySource.parent_hashes:semantic,SourceColorWaiverSubjectIdentitySource.expires_day:semantic,SourceColorWaiverSubjectIdentitySource.signature:nonsemantic:self-signature-is-outside-its-subject",
    "source_bindings=SourceColorWaiverSubjectIdentitySource.key_id>key-id-byte-count+key-id,SourceColorWaiverSubjectIdentitySource.scope>scope-byte-count+scope,SourceColorWaiverSubjectIdentitySource.node_name>node-name-byte-count+node-name,SourceColorWaiverSubjectIdentitySource.claimed_color>claimed-color-byte-count+claimed-color-canonical-bytes,SourceColorWaiverSubjectIdentitySource.annotation_id>annotation-id-byte-count+annotation-id,SourceColorWaiverSubjectIdentitySource.annotation_signer>annotation-signer-byte-count+annotation-signer,SourceColorWaiverSubjectIdentitySource.annotation_reason>annotation-reason-byte-count+annotation-reason,SourceColorWaiverSubjectIdentitySource.parent_hashes>parent-count+parent-order+parent-hashes,SourceColorWaiverSubjectIdentitySource.expires_day>expires-day",
    "external_semantic_fields=transport-version,domain-byte-count,preimage-domain,source-operation-sentinel",
    "semantic_fields=transport-version,domain-byte-count,preimage-domain,source-operation-sentinel,key-id-byte-count,key-id,scope-byte-count,scope,node-name-byte-count,node-name,claimed-color-byte-count,claimed-color-canonical-bytes,annotation-id-byte-count,annotation-id,annotation-signer-byte-count,annotation-signer,annotation-reason-byte-count,annotation-reason,parent-count,parent-order,parent-hashes,expires-day",
    "excluded_fields=admission-day:verification-context-not-signed-subject,policy-fingerprint:verifier-result-not-signed-subject",
    "consumers=WaiverGrant::signing_payload_source,WaiverVerifier::verify,ColorGraph::source_waived,ColorGraph::verify_replay",
    "mutations=transport-version:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_identity_fields_move_independently,domain-byte-count:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_identity_fields_move_independently,preimage-domain:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_identity_fields_move_independently,source-operation-sentinel:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_identity_fields_move_independently,key-id-byte-count:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_identity_fields_move_independently,key-id:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_identity_fields_move_independently,scope-byte-count:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_identity_fields_move_independently,scope:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_identity_fields_move_independently,node-name-byte-count:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_identity_fields_move_independently,node-name:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_identity_fields_move_independently,claimed-color-byte-count:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_identity_fields_move_independently,claimed-color-canonical-bytes:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_identity_fields_move_independently,annotation-id-byte-count:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_identity_fields_move_independently,annotation-id:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_identity_fields_move_independently,annotation-signer-byte-count:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_identity_fields_move_independently,annotation-signer:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_identity_fields_move_independently,annotation-reason-byte-count:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_identity_fields_move_independently,annotation-reason:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_identity_fields_move_independently,parent-count:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_identity_fields_move_independently,parent-order:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_identity_fields_move_independently,parent-hashes:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_identity_fields_move_independently,expires-day:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_identity_fields_move_independently",
    "nonsemantic_mutations=SourceColorWaiverSubjectIdentitySource.signature:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_excluded_fields_do_not_move_identity,admission-day:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_excluded_fields_do_not_move_identity,policy-fingerprint:crates/fs-ledger/src/lib.rs#source_color_waiver_subject_excluded_fields_do_not_move_identity",
    "field_guard=classify_source_color_waiver_subject_identity_fields",
    "transport_guard=ledger_source_color_waiver_subject_transport_guard",
    "version_guard=crates/fs-ledger/src/lib.rs#ledger_semantic_identity_versions_fail_closed",
    "coupling_surface=fs-ledger:source-color-waiver-subject",
];

/// Owner declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const COLOR_NODE_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-ledger:color-node",
    "version_const=COLOR_NODE_IDENTITY_VERSION",
    "version=9",
    "domain=org.frankensim.fs-ledger.color-node.v9",
    "domain_const=COLOR_NODE_IDENTITY_DOMAIN",
    "encoder=ledger_color_node_identity",
    "encoder_helpers=ledger_color_node_identity_with_schema",
    "schema_constants=COLOR_NODE_IDENTITY_VERSION,COLOR_NODE_IDENTITY_DOMAIN,COLOR_NODE_PREIMAGE_DOMAIN,crates/fs-ledger/src/colors.rs#COLOR_NODE_HASH_ENCODING_VERSION,crates/fs-ledger/src/colors.rs#COLOR_NODE_HASH_DOMAIN,crates/fs-ledger/src/colors.rs#COLOR_WRITE_ROW_SCHEMA_VERSION,crates/fs-evidence/src/color.rs#COLOR_ALGEBRA_VERSION",
    "schema_functions=crates/fs-ledger/src/colors.rs#ColorGraph::node_hash,crates/fs-ledger/src/colors.rs#ColorGraph::node_hash_from_canonical_payloads,crates/fs-ledger/src/colors.rs#source_origin_canonical_bytes,crates/fs-ledger/src/colors.rs#WaiverGrant::signing_payload_for,crates/fs-ledger/src/colors.rs#push_len,crates/fs-ledger/src/colors.rs#push_field,crates/fs-ledger/src/colors.rs#interval_op_tag,crates/fs-evidence/src/color.rs#Color::canonical_bytes,identity_schema_is_current",
    "schema_dependencies=fs-ledger:derived-color-waiver-subject,fs-ledger:source-color-waiver-subject,fs-ledger:source-origin-request",
    "digest=blake3-256-over-canonical-provenance-transport",
    "encoding=typed-binary",
    "sources=ColorNodeIdentitySource",
    "source_fields=ColorNodeIdentitySource.node_id:nonsemantic:ledger-local-row-identity,ColorNodeIdentitySource.operation_tag:semantic,ColorNodeIdentitySource.name:semantic,ColorNodeIdentitySource.color:semantic,ColorNodeIdentitySource.parent_local_ids:nonsemantic:ledger-local-parent-addresses,ColorNodeIdentitySource.parent_hashes:semantic,ColorNodeIdentitySource.demotions:semantic,ColorNodeIdentitySource.origin:semantic,ColorNodeIdentitySource.origin_policy_fingerprint:semantic,ColorNodeIdentitySource.waiver_dependencies:semantic,ColorNodeIdentitySource.waiver:semantic,ColorNodeIdentitySource.grant_payload:semantic,ColorNodeIdentitySource.grant_signature:semantic,ColorNodeIdentitySource.waiver_policy_fingerprint:semantic,ColorNodeIdentitySource.waiver_admission_day:semantic,ColorNodeIdentitySource.stored_hash:derived:recomputed-identity-output",
    "source_bindings=ColorNodeIdentitySource.operation_tag>operation-presence+operation-tag,ColorNodeIdentitySource.name>node-name-byte-count+node-name,ColorNodeIdentitySource.color>color-byte-count+color-canonical-bytes,ColorNodeIdentitySource.parent_hashes>parent-count+parent-order+parent-hashes,ColorNodeIdentitySource.demotions>demotion-count+demotion-order+demotion-canonical-bytes,ColorNodeIdentitySource.origin>origin-presence+origin-canonical-bytes,ColorNodeIdentitySource.origin_policy_fingerprint>origin-policy-presence+origin-policy-fingerprint,ColorNodeIdentitySource.waiver_dependencies>waiver-dependency-count+waiver-dependency-order+waiver-dependency-canonical-bytes,ColorNodeIdentitySource.waiver>waiver-presence+waiver-canonical-bytes,ColorNodeIdentitySource.grant_payload>grant-presence+grant-payload,ColorNodeIdentitySource.grant_signature>grant-signature,ColorNodeIdentitySource.waiver_policy_fingerprint>waiver-policy-presence+waiver-policy-fingerprint,ColorNodeIdentitySource.waiver_admission_day>waiver-admission-day-presence+waiver-admission-day",
    "external_semantic_fields=transport-version,domain-byte-count,preimage-domain",
    "semantic_fields=transport-version,domain-byte-count,preimage-domain,operation-presence,operation-tag,node-name-byte-count,node-name,color-byte-count,color-canonical-bytes,parent-count,parent-order,parent-hashes,demotion-count,demotion-order,demotion-canonical-bytes,origin-presence,origin-canonical-bytes,origin-policy-presence,origin-policy-fingerprint,waiver-dependency-count,waiver-dependency-order,waiver-dependency-canonical-bytes,waiver-presence,waiver-canonical-bytes,grant-presence,grant-payload,grant-signature,waiver-policy-presence,waiver-policy-fingerprint,waiver-admission-day-presence,waiver-admission-day",
    "excluded_fields=display-json:audit-rendering-is-not-canonical,write-timestamp:wall-clock-envelope,color-row-schema:storage-envelope-version",
    "consumers=ColorGraph::source,ColorGraph::source_with_origin,ColorGraph::source_waived,ColorGraph::derive,ColorGraph::derive_waived,ColorGraph::verify_replay,ColorNode::hash",
    "mutations=transport-version:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,domain-byte-count:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,preimage-domain:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,operation-presence:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,operation-tag:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,node-name-byte-count:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,node-name:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,color-byte-count:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,color-canonical-bytes:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,parent-count:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,parent-order:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,parent-hashes:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,demotion-count:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,demotion-order:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,demotion-canonical-bytes:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,origin-presence:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,origin-canonical-bytes:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,origin-policy-presence:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,origin-policy-fingerprint:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,waiver-dependency-count:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,waiver-dependency-order:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,waiver-dependency-canonical-bytes:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,waiver-presence:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,waiver-canonical-bytes:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,grant-presence:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,grant-payload:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,grant-signature:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,waiver-policy-presence:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,waiver-policy-fingerprint:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,waiver-admission-day-presence:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently,waiver-admission-day:crates/fs-ledger/src/lib.rs#color_node_identity_fields_move_independently",
    "nonsemantic_mutations=ColorNodeIdentitySource.node_id:crates/fs-ledger/src/lib.rs#color_node_excluded_fields_do_not_move_identity,ColorNodeIdentitySource.parent_local_ids:crates/fs-ledger/src/lib.rs#color_node_excluded_fields_do_not_move_identity,display-json:crates/fs-ledger/src/lib.rs#color_node_excluded_fields_do_not_move_identity,write-timestamp:crates/fs-ledger/src/lib.rs#color_node_excluded_fields_do_not_move_identity,color-row-schema:crates/fs-ledger/src/lib.rs#color_node_excluded_fields_do_not_move_identity",
    "field_guard=classify_color_node_identity_fields",
    "transport_guard=ledger_color_node_identity",
    "version_guard=crates/fs-ledger/src/lib.rs#ledger_semantic_identity_versions_fail_closed",
    "coupling_surface=fs-ledger:color-node",
];

/// Owner declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const COLOR_ADMISSION_POLICY_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-ledger:color-admission-policy",
    "version_const=COLOR_ADMISSION_POLICY_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-ledger.color-admission-policy.v1",
    "domain_const=COLOR_ADMISSION_POLICY_IDENTITY_DOMAIN",
    "encoder=ledger_color_admission_policy_identity",
    "encoder_helpers=ledger_color_admission_policy_identity_with_schema",
    "schema_constants=COLOR_ADMISSION_POLICY_IDENTITY_VERSION,COLOR_ADMISSION_POLICY_IDENTITY_DOMAIN,COLOR_ADMISSION_POLICY_PREIMAGE_DOMAIN,crates/fs-ledger/src/colors.rs#COLOR_WRITE_ROW_SCHEMA_VERSION,crates/fs-evidence/src/color.rs#COLOR_ALGEBRA_VERSION",
    "schema_functions=crates/fs-ledger/src/colors.rs#color_admission_policy_fingerprint,crates/fs-ledger/src/colors.rs#ColorGraph::admission_receipt,crates/fs-ledger/src/colors.rs#ColorGraph::admission_receipt_in_regime,crates/fs-ledger/src/colors.rs#LedgerColorAdmissionVerifier::verify,identity_schema_is_current",
    "schema_dependencies=fs-ledger:color-node",
    "digest=blake3-256-policy-fingerprint",
    "encoding=typed-binary",
    "sources=ColorAdmissionPolicyIdentitySource",
    "source_fields=ColorAdmissionPolicyIdentitySource.color_write_row_schema_version:semantic,ColorAdmissionPolicyIdentitySource.color_algebra_version:semantic",
    "source_bindings=ColorAdmissionPolicyIdentitySource.color_write_row_schema_version>color-write-row-schema-version,ColorAdmissionPolicyIdentitySource.color_algebra_version>color-algebra-version",
    "external_semantic_fields=preimage-domain",
    "semantic_fields=preimage-domain,color-write-row-schema-version,color-algebra-version",
    "excluded_fields=build-version:build-envelope-not-policy-semantics,wall-clock:admission-time-not-policy-identity",
    "consumers=color_admission_policy_fingerprint,ColorGraph::admission_receipt,LedgerColorAdmissionVerifier::verify,fs-evidence::AdmittedColor",
    "mutations=preimage-domain:crates/fs-ledger/src/lib.rs#color_admission_policy_identity_fields_move_independently,color-write-row-schema-version:crates/fs-ledger/src/lib.rs#color_admission_policy_identity_fields_move_independently,color-algebra-version:crates/fs-ledger/src/lib.rs#color_admission_policy_identity_fields_move_independently",
    "nonsemantic_mutations=build-version:crates/fs-ledger/src/lib.rs#color_admission_policy_excluded_fields_do_not_move_identity,wall-clock:crates/fs-ledger/src/lib.rs#color_admission_policy_excluded_fields_do_not_move_identity",
    "field_guard=classify_color_admission_policy_identity_fields",
    "transport_guard=ledger_color_admission_policy_identity",
    "version_guard=crates/fs-ledger/src/lib.rs#ledger_semantic_identity_versions_fail_closed",
    "coupling_surface=fs-ledger:color-admission-policy",
];

/// Owner declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const VCS_LEDGER_LINEAGE_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-ledger:vcs-ledger-lineage",
    "version_const=VCS_LEDGER_LINEAGE_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-ledger.vcs-ledger-lineage.v1",
    "domain_const=VCS_LEDGER_LINEAGE_IDENTITY_DOMAIN",
    "encoder=ledger_vcs_ledger_lineage_identity",
    "encoder_helpers=ledger_vcs_ledger_lineage_identity_with_domain",
    "schema_constants=VCS_LEDGER_LINEAGE_IDENTITY_VERSION,VCS_LEDGER_LINEAGE_IDENTITY_DOMAIN,VCS_LEDGER_LINEAGE_PREIMAGE_DOMAIN,crates/fs-ledger/src/vcs.rs#LEDGER_IDENTITY_DOMAIN,VCS_IDENTITY_EVENT_KIND",
    "schema_functions=crates/fs-ledger/src/vcs.rs#hash_frame,crates/fs-ledger/src/vcs.rs#domain_hasher,crates/fs-ledger/src/vcs.rs#hash_field,crates/fs-ledger/src/vcs.rs#framed_hash,crates/fs-ledger/src/vcs.rs#Ledger::vcs_identity,Ledger::append_vcs_identity_event,identity_schema_is_current",
    "schema_dependencies=none",
    "digest=blake3-256-domain-separated",
    "encoding=typed-binary",
    "sources=VcsLedgerLineageIdentitySource",
    "source_fields=VcsLedgerLineageIdentitySource.mint_path:semantic,VcsLedgerLineageIdentitySource.minted_ns:semantic",
    "source_bindings=VcsLedgerLineageIdentitySource.mint_path>mint-path-byte-count+mint-path,VcsLedgerLineageIdentitySource.minted_ns>minted-nanoseconds",
    "external_semantic_fields=domain-label-frame,domain-byte-count,preimage-domain",
    "semantic_fields=domain-label-frame,domain-byte-count,preimage-domain,mint-path-byte-count,mint-path,minted-nanoseconds",
    "excluded_fields=identity-event-rowid:storage-envelope-only,current-path-after-mint:persisted-lineage-survives-moves,reopen-time:reopen-is-not-remint",
    "consumers=Ledger::vcs_identity,Vcs::commit,Vcs::lookup,Vcs::checkout,fs-ledger:vcs-commit-envelope",
    "mutations=domain-label-frame:crates/fs-ledger/src/lib.rs#vcs_ledger_lineage_identity_fields_move_independently,domain-byte-count:crates/fs-ledger/src/lib.rs#vcs_ledger_lineage_identity_fields_move_independently,preimage-domain:crates/fs-ledger/src/lib.rs#vcs_ledger_lineage_identity_fields_move_independently,mint-path-byte-count:crates/fs-ledger/src/lib.rs#vcs_ledger_lineage_identity_fields_move_independently,mint-path:crates/fs-ledger/src/lib.rs#vcs_ledger_lineage_identity_fields_move_independently,minted-nanoseconds:crates/fs-ledger/src/lib.rs#vcs_ledger_lineage_identity_fields_move_independently",
    "nonsemantic_mutations=identity-event-rowid:crates/fs-ledger/src/lib.rs#vcs_ledger_lineage_excluded_fields_do_not_move_identity,current-path-after-mint:crates/fs-ledger/src/lib.rs#vcs_ledger_lineage_excluded_fields_do_not_move_identity,reopen-time:crates/fs-ledger/src/lib.rs#vcs_ledger_lineage_excluded_fields_do_not_move_identity",
    "field_guard=classify_vcs_ledger_lineage_identity_fields",
    "transport_guard=ledger_vcs_ledger_lineage_identity",
    "version_guard=crates/fs-ledger/src/lib.rs#ledger_semantic_identity_versions_fail_closed",
    "coupling_surface=fs-ledger:vcs-ledger-lineage",
];

/// Owner declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const VCS_COMMIT_LEAF_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-ledger:vcs-commit-leaf",
    "version_const=VCS_COMMIT_LEAF_IDENTITY_VERSION",
    "version=2",
    "domain=org.frankensim.fs-ledger.vcs-commit-leaf.v2",
    "domain_const=VCS_COMMIT_LEAF_IDENTITY_DOMAIN",
    "encoder=ledger_vcs_commit_leaf_identity",
    "encoder_helpers=ledger_vcs_commit_leaf_identity_with_domain,identity_vcs_hash_optional_field",
    "schema_constants=VCS_COMMIT_LEAF_IDENTITY_VERSION,VCS_COMMIT_LEAF_IDENTITY_DOMAIN,VCS_COMMIT_LEAF_PREIMAGE_DOMAIN,crates/fs-ledger/src/vcs.rs#COMMIT_LEAF_DOMAIN",
    "schema_functions=crates/fs-ledger/src/vcs.rs#hash_frame,crates/fs-ledger/src/vcs.rs#domain_hasher,crates/fs-ledger/src/vcs.rs#hash_field,crates/fs-ledger/src/vcs.rs#hash_optional_field,crates/fs-ledger/src/vcs.rs#Ledger::commit_leaf,crates/fs-ledger/src/vcs.rs#Ledger::op_artifact_edges,crates/fs-ledger/src/vcs.rs#Ledger::commit_exec_mode,identity_schema_is_current",
    "schema_dependencies=fs-ledger:artifact-content",
    "digest=blake3-256-domain-separated",
    "encoding=typed-binary",
    "sources=VcsCommitLeafIdentitySource,VcsCommitEdgeIdentitySource",
    "source_fields=VcsCommitLeafIdentitySource.ir:semantic,VcsCommitLeafIdentitySource.seed:semantic,VcsCommitLeafIdentitySource.versions:semantic,VcsCommitLeafIdentitySource.budget:semantic,VcsCommitLeafIdentitySource.capability:semantic,VcsCommitLeafIdentitySource.outcome:semantic,VcsCommitLeafIdentitySource.diagnostic:semantic,VcsCommitLeafIdentitySource.execution_mode:semantic,VcsCommitLeafIdentitySource.edges:semantic,VcsCommitEdgeIdentitySource.role:semantic,VcsCommitEdgeIdentitySource.artifact_hash:semantic",
    "source_bindings=VcsCommitLeafIdentitySource.ir>ir-byte-count+ir-bytes,VcsCommitLeafIdentitySource.seed>seed-byte-count+seed-bytes,VcsCommitLeafIdentitySource.versions>versions-byte-count+versions-bytes,VcsCommitLeafIdentitySource.budget>budget-byte-count+budget-bytes,VcsCommitLeafIdentitySource.capability>capability-byte-count+capability-bytes,VcsCommitLeafIdentitySource.outcome>outcome-presence+outcome-byte-count+outcome-bytes,VcsCommitLeafIdentitySource.diagnostic>diagnostic-presence+diagnostic-byte-count+diagnostic-bytes,VcsCommitLeafIdentitySource.execution_mode>execution-mode-byte-count+execution-mode,VcsCommitLeafIdentitySource.edges>edge-count+edge-order,VcsCommitEdgeIdentitySource.role>edge-role-byte-count+edge-role,VcsCommitEdgeIdentitySource.artifact_hash>artifact-hash",
    "external_semantic_fields=domain-label-frame,domain-byte-count,preimage-domain",
    "semantic_fields=domain-label-frame,domain-byte-count,preimage-domain,ir-byte-count,ir-bytes,seed-byte-count,seed-bytes,versions-byte-count,versions-bytes,budget-byte-count,budget-bytes,capability-byte-count,capability-bytes,outcome-presence,outcome-byte-count,outcome-bytes,diagnostic-presence,diagnostic-byte-count,diagnostic-bytes,execution-mode-byte-count,execution-mode,edge-count,edge-order,edge-role-byte-count,edge-role,artifact-hash",
    "excluded_fields=op-rowid:ledger-local-address,session:execution-envelope,t-start:wall-clock-envelope,t-end:wall-clock-envelope,branch-id:ledger-local-envelope,edge-rowid:storage-envelope-only",
    "consumers=Ledger::commit_leaf,Vcs::commit,Vcs::checkout_delta,Vcs::lookup_semantic,fs-ledger:vcs-commit-root",
    "mutations=domain-label-frame:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,domain-byte-count:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,preimage-domain:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,ir-byte-count:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,ir-bytes:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,seed-byte-count:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,seed-bytes:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,versions-byte-count:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,versions-bytes:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,budget-byte-count:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,budget-bytes:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,capability-byte-count:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,capability-bytes:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,outcome-presence:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,outcome-byte-count:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,outcome-bytes:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,diagnostic-presence:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,diagnostic-byte-count:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,diagnostic-bytes:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,execution-mode-byte-count:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,execution-mode:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,edge-count:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,edge-order:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,edge-role-byte-count:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,edge-role:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently,artifact-hash:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_identity_fields_move_independently",
    "nonsemantic_mutations=op-rowid:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_excluded_fields_do_not_move_identity,session:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_excluded_fields_do_not_move_identity,t-start:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_excluded_fields_do_not_move_identity,t-end:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_excluded_fields_do_not_move_identity,branch-id:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_excluded_fields_do_not_move_identity,edge-rowid:crates/fs-ledger/src/lib.rs#vcs_commit_leaf_excluded_fields_do_not_move_identity",
    "field_guard=classify_vcs_commit_leaf_identity_fields",
    "transport_guard=ledger_vcs_commit_leaf_identity",
    "version_guard=crates/fs-ledger/src/lib.rs#ledger_semantic_identity_versions_fail_closed",
    "coupling_surface=fs-ledger:vcs-commit-leaf",
];

/// Owner declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const VCS_COMMIT_ROOT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-ledger:vcs-commit-root",
    "version_const=VCS_COMMIT_ROOT_IDENTITY_VERSION",
    "version=2",
    "domain=org.frankensim.fs-ledger.vcs-commit-root.v2",
    "domain_const=VCS_COMMIT_ROOT_IDENTITY_DOMAIN",
    "encoder=ledger_vcs_commit_root_identity",
    "encoder_helpers=ledger_vcs_commit_root_identity_with_domains,identity_vcs_hash_frame,identity_vcs_domain_hasher,identity_vcs_hash_field,identity_vcs_framed_hash",
    "schema_constants=VCS_COMMIT_ROOT_IDENTITY_VERSION,VCS_COMMIT_ROOT_IDENTITY_DOMAIN,VCS_MERKLE_PAIR_PREIMAGE_DOMAIN,VCS_MERKLE_ODD_PREIMAGE_DOMAIN,VCS_COMMIT_ROOT_PREIMAGE_DOMAIN,crates/fs-ledger/src/vcs.rs#MERKLE_PAIR_DOMAIN,crates/fs-ledger/src/vcs.rs#MERKLE_ODD_DOMAIN,crates/fs-ledger/src/vcs.rs#COMMIT_ROOT_DOMAIN",
    "schema_functions=crates/fs-ledger/src/vcs.rs#hash_frame,crates/fs-ledger/src/vcs.rs#domain_hasher,crates/fs-ledger/src/vcs.rs#hash_field,crates/fs-ledger/src/vcs.rs#framed_hash,crates/fs-ledger/src/vcs.rs#merkle_root,crates/fs-ledger/src/vcs.rs#Vcs::commit,identity_schema_is_current",
    "schema_dependencies=fs-ledger:vcs-commit-leaf",
    "digest=blake3-256-domain-separated-binary-merkle",
    "encoding=typed-binary",
    "sources=VcsCommitRootIdentitySource",
    "source_fields=VcsCommitRootIdentitySource.leaves:semantic",
    "source_bindings=VcsCommitRootIdentitySource.leaves>leaf-count+leaf-order+leaf-hashes+tree-shape",
    "external_semantic_fields=merkle-domain-set",
    "semantic_fields=merkle-domain-set,leaf-count,leaf-order,leaf-hashes,tree-shape",
    "excluded_fields=ledger-identity:semantic-state-is-portable,branch-id:semantic-state-is-portable,local-op-ids:leaf-hashes-are-portable,commit-time:wall-clock-envelope",
    "consumers=Vcs::commit,Vcs::lookup_semantic,Vcs::checkout_delta,fs-ledger:vcs-commit-envelope",
    "mutations=merkle-domain-set:crates/fs-ledger/src/lib.rs#vcs_commit_root_identity_fields_move_independently,leaf-count:crates/fs-ledger/src/lib.rs#vcs_commit_root_identity_fields_move_independently,leaf-order:crates/fs-ledger/src/lib.rs#vcs_commit_root_identity_fields_move_independently,leaf-hashes:crates/fs-ledger/src/lib.rs#vcs_commit_root_identity_fields_move_independently,tree-shape:crates/fs-ledger/src/lib.rs#vcs_commit_root_identity_fields_move_independently",
    "nonsemantic_mutations=ledger-identity:crates/fs-ledger/src/lib.rs#vcs_commit_root_excluded_fields_do_not_move_identity,branch-id:crates/fs-ledger/src/lib.rs#vcs_commit_root_excluded_fields_do_not_move_identity,local-op-ids:crates/fs-ledger/src/lib.rs#vcs_commit_root_excluded_fields_do_not_move_identity,commit-time:crates/fs-ledger/src/lib.rs#vcs_commit_root_excluded_fields_do_not_move_identity",
    "field_guard=classify_vcs_commit_root_identity_fields",
    "transport_guard=ledger_vcs_commit_root_identity",
    "version_guard=crates/fs-ledger/src/lib.rs#ledger_semantic_identity_versions_fail_closed",
    "coupling_surface=fs-ledger:vcs-commit-root",
];

/// Owner declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const VCS_COMMIT_ENVELOPE_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-ledger:vcs-commit-envelope",
    "version_const=VCS_COMMIT_ENVELOPE_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-ledger.vcs-commit-envelope.v1",
    "domain_const=VCS_COMMIT_ENVELOPE_IDENTITY_DOMAIN",
    "encoder=ledger_vcs_commit_envelope_identity",
    "encoder_helpers=none",
    "schema_constants=VCS_COMMIT_ENVELOPE_IDENTITY_VERSION,VCS_COMMIT_ENVELOPE_IDENTITY_DOMAIN",
    "schema_functions=crates/fs-ledger/src/vcs.rs#CommitInfo::id,crates/fs-ledger/src/vcs.rs#Vcs::commit,crates/fs-ledger/src/vcs.rs#Vcs::lookup,crates/fs-ledger/src/vcs.rs#Vcs::checkout,identity_schema_is_current",
    "schema_dependencies=fs-ledger:vcs-commit-root,fs-ledger:vcs-ledger-lineage",
    "digest=none-fixed-width-envelope-key",
    "encoding=fixed-width-key",
    "sources=VcsCommitEnvelopeIdentitySource",
    "source_fields=VcsCommitEnvelopeIdentitySource.ledger:semantic,VcsCommitEnvelopeIdentitySource.branch:semantic,VcsCommitEnvelopeIdentitySource.root:semantic",
    "source_bindings=VcsCommitEnvelopeIdentitySource.ledger>ledger-lineage,VcsCommitEnvelopeIdentitySource.branch>branch-id,VcsCommitEnvelopeIdentitySource.root>semantic-root",
    "external_semantic_fields=none",
    "semantic_fields=ledger-lineage,branch-id,semantic-root",
    "excluded_fields=frontier-op:ledger-local-snapshot-envelope,parent-root:commit-graph-edge-not-key,event-rowid:storage-envelope-only,commit-time:wall-clock-envelope",
    "consumers=CommitInfo::id,Vcs::commit,Vcs::lookup,Vcs::checkout,Vcs::checkout_delta",
    "mutations=ledger-lineage:crates/fs-ledger/src/lib.rs#vcs_commit_envelope_identity_fields_move_independently,branch-id:crates/fs-ledger/src/lib.rs#vcs_commit_envelope_identity_fields_move_independently,semantic-root:crates/fs-ledger/src/lib.rs#vcs_commit_envelope_identity_fields_move_independently",
    "nonsemantic_mutations=frontier-op:crates/fs-ledger/src/lib.rs#vcs_commit_envelope_excluded_fields_do_not_move_identity,parent-root:crates/fs-ledger/src/lib.rs#vcs_commit_envelope_excluded_fields_do_not_move_identity,event-rowid:crates/fs-ledger/src/lib.rs#vcs_commit_envelope_excluded_fields_do_not_move_identity,commit-time:crates/fs-ledger/src/lib.rs#vcs_commit_envelope_excluded_fields_do_not_move_identity",
    "field_guard=classify_vcs_commit_envelope_identity_fields",
    "transport_guard=ledger_vcs_commit_envelope_identity",
    "version_guard=crates/fs-ledger/src/lib.rs#ledger_semantic_identity_versions_fail_closed",
    "coupling_surface=fs-ledger:vcs-commit-envelope",
];

pub(crate) const VCS_IDENTITY_EVENT_KIND: &str = "vcs-identity";

/// Opaque identity of one physical ledger instance.
///
/// File-backed ledgers persist this value in schema metadata, so aliases and
/// reopenings agree while a replacement database at the same path does not.
/// In-memory ledgers retain a generated value inside the handle, so moving the
/// Rust value cannot change its identity and independent handles never alias by
/// address reuse.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LedgerInstanceId([u8; 16]);

impl LedgerInstanceId {
    /// Exact RFC 4122-shaped UUID bytes carrying the opaque identity.
    #[must_use]
    pub const fn as_bytes(self) -> [u8; 16] {
        self.0
    }

    /// Stable lowercase UUID rendering for diagnostics and manifests.
    #[must_use]
    pub fn to_uuid_string(self) -> String {
        self.to_string()
    }
}

impl std::fmt::Display for LedgerInstanceId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let bytes = self.0;
        write!(
            f,
            "{:02x}{:02x}{:02x}{:02x}-{:02x}{:02x}-{:02x}{:02x}-{:02x}{:02x}-\
             {:02x}{:02x}{:02x}{:02x}{:02x}{:02x}",
            bytes[0],
            bytes[1],
            bytes[2],
            bytes[3],
            bytes[4],
            bytes[5],
            bytes[6],
            bytes[7],
            bytes[8],
            bytes[9],
            bytes[10],
            bytes[11],
            bytes[12],
            bytes[13],
            bytes[14],
            bytes[15],
        )
    }
}

impl std::fmt::Debug for LedgerInstanceId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "LedgerInstanceId({self})")
    }
}

fn fresh_ledger_instance_id() -> Result<LedgerInstanceId, LedgerError> {
    #[cfg(unix)]
    {
        use std::io::Read as _;

        let mut uuid = [0_u8; 16];
        let mut entropy = std::fs::File::open("/dev/urandom").map_err(|error| {
            LedgerError::InstanceIdentityUnavailable {
                detail: format!("cannot open /dev/urandom: {error}"),
            }
        })?;
        entropy.read_exact(&mut uuid).map_err(|error| {
            LedgerError::InstanceIdentityUnavailable {
                detail: format!("cannot read 16 bytes from /dev/urandom: {error}"),
            }
        })?;
        uuid[6] = (uuid[6] & 0x0f) | 0x40;
        uuid[8] = (uuid[8] & 0x3f) | 0x80;
        Ok(LedgerInstanceId(uuid))
    }

    #[cfg(not(unix))]
    {
        Err(LedgerError::InstanceIdentityUnavailable {
            detail: format!(
                "no safe std-only OS entropy source is implemented for target OS {}",
                std::env::consts::OS
            ),
        })
    }
}

fn decode_ledger_instance_id(
    row_count: usize,
    singleton: Option<&SqliteValue>,
    instance_id: Option<&SqliteValue>,
) -> Result<LedgerInstanceId, LedgerError> {
    if row_count != 1 {
        return Err(LedgerError::InstanceIdentityCorrupt {
            detail: format!(
                "ledger_identity must contain exactly one total row, found {row_count}"
            ),
        });
    }
    if !matches!(singleton, Some(SqliteValue::Integer(1))) {
        return Err(LedgerError::InstanceIdentityCorrupt {
            detail: format!("ledger_identity singleton key must be INTEGER 1, found {singleton:?}"),
        });
    }
    let Some(SqliteValue::Blob(bytes)) = instance_id else {
        return Err(LedgerError::InstanceIdentityCorrupt {
            detail: "ledger_identity.instance_id is not a BLOB".to_string(),
        });
    };
    let uuid: [u8; 16] =
        bytes
            .as_ref()
            .try_into()
            .map_err(|_| LedgerError::InstanceIdentityCorrupt {
                detail: format!(
                    "ledger_identity.instance_id must be a 16-byte UUID, found {} bytes",
                    bytes.len()
                ),
            })?;
    if uuid[6] & 0xf0 != 0x40 || uuid[8] & 0xc0 != 0x80 {
        return Err(LedgerError::InstanceIdentityCorrupt {
            detail: "ledger_identity.instance_id has invalid UUID version or variant bits"
                .to_string(),
        });
    }
    Ok(LedgerInstanceId(uuid))
}

// ---------------------------------------------------------------------------
// Error model
// ---------------------------------------------------------------------------

/// Structured, machine-actionable ledger errors (Decalogue P10). Never a
/// panic across the crate boundary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LedgerError {
    /// The database could not be opened or configured.
    Open {
        /// Path passed to [`Ledger::open`].
        path: String,
        /// Underlying engine detail.
        detail: String,
    },
    /// The file was written by a newer fs-ledger than this build supports.
    FutureSchema {
        /// `PRAGMA user_version` found in the file.
        found: i64,
        /// Highest version this build can read/write.
        supported: i64,
    },
    /// A SQL-level failure that is not retryable contention.
    Sql {
        /// What the ledger was doing.
        context: String,
        /// Underlying engine detail.
        detail: String,
    },
    /// Retryable contention (busy/locked/write-conflict). Retry the call,
    /// ideally with backoff; no state was silently lost.
    Busy {
        /// What the ledger was doing.
        context: String,
        /// Underlying engine detail.
        detail: String,
    },
    /// A required Five Explicits field is missing or malformed (P4/P10).
    MissingExplicit {
        /// Field name: `seed`, `versions`, `budget`, or `capability`.
        field: String,
        /// What is wrong and how to fix it.
        problem: String,
    },
    /// An input failed validation (structured; names the field).
    Invalid {
        /// Field name.
        field: String,
        /// What is wrong and how to fix it.
        problem: String,
    },
    /// A retained state checkpoint names law semantics that differ from the
    /// exact tuple supplied by the replayer. State bytes are withheld.
    UnknownStateSemantics {
        /// Durable state-slot digest rendered as lowercase hex.
        state_slot_hex: String,
        /// Law id stored with the checkpoint.
        stored_law: String,
        /// Stored law version.
        stored_law_version: u32,
        /// Stored state-schema version.
        stored_state_schema_version: u32,
        /// Law id the caller knows how to replay.
        expected_law: String,
        /// Expected law version.
        expected_law_version: u32,
        /// Expected state-schema version.
        expected_state_schema_version: u32,
        /// Stored L1 canonical parameter-block hash.
        stored_parameters_hash: String,
        /// Expected L1 canonical parameter-block hash.
        expected_parameters_hash: String,
        /// Stored L3 contract plus implementation hash.
        stored_contract_and_code_hash: String,
        /// Expected L3 contract plus implementation hash.
        expected_contract_and_code_hash: String,
        /// Closed semantic fields that disagreed.
        differences: Vec<&'static str>,
    },
    /// Stored bytes no longer match their recorded content hash.
    Corrupt {
        /// Hex hash of the first corrupted artifact.
        hash_hex: String,
        /// Diagnosis.
        detail: String,
    },
    /// Stored metadata declares an artifact larger than the caller's explicit
    /// materialization budget. The payload is refused before any byte callback
    /// or allocation; this refusal makes no independent integrity claim.
    ArtifactReadLimit {
        /// Hex content hash of the refused artifact.
        hash_hex: String,
        /// Caller-supplied maximum payload bytes.
        limit: u64,
        /// Stored payload length observed during metadata-only preflight.
        observed: u64,
    },
    /// An operation row bypassed the canonical write path and violates the
    /// bounded storage contract. Reads refuse it before materializing any
    /// variable-size field.
    OpCorrupt {
        /// Operation id addressed by the read.
        op: i64,
        /// Bounded structural diagnosis without hostile stored values.
        detail: String,
    },
    /// A tune row bypassed the canonical write path and violates the bounded
    /// storage contract. Reads refuse the row rather than materializing or
    /// interpreting it.
    TuneCorrupt {
        /// Canonical kernel identity used to address the row or scan.
        kernel: String,
        /// Diagnosis without embedding hostile stored values.
        detail: String,
    },
    /// A tune history exceeds a deterministic read budget. The caller must
    /// narrow or compact the history before retrying.
    TuneReadLimit {
        /// Canonical kernel identity whose history was refused.
        kernel: String,
        /// Bounded resource: `rows` or `materialized_bytes`.
        resource: &'static str,
        /// Configured maximum.
        limit: usize,
        /// Exact observation or a conservative lower bound.
        observed_at_least: usize,
    },
    /// A referenced row does not exist.
    NotFound {
        /// Description of the missing row.
        what: String,
    },
    /// `finish_op` was called on an op that already has an outcome.
    DoubleFinish {
        /// The op id.
        op: i64,
    },
    /// A lineage edge was offered after the operation reached its terminal
    /// outcome. `finish_op` is the public-API sealing transition: callers must
    /// attach every input and output before publishing the outcome.
    OpLineageSealed {
        /// The terminal operation whose edge set cannot be extended.
        op: i64,
    },
    /// An [`ArtifactWriter`] cannot be opened inside an explicit transaction
    /// (it owns its own transaction for crash atomicity).
    WriterInTransaction,
    /// The on-disk objects do not match the schema its `user_version`
    /// claims (bead gp3.18): pre-existing incompatible tables, a
    /// partially initialized file, wrong columns/affinities, missing
    /// indexes, foreign objects, or an altered identity guard. Refused
    /// before migration advances the version or a checked authority boundary
    /// trusts the identity — `CREATE TABLE IF NOT EXISTS` must never launder
    /// an alien schema into a labeled one.
    SchemaMismatch {
        /// The `PRAGMA user_version` the file claims.
        claimed_version: i64,
        /// Every attestation violation found (object-level diffs).
        violations: Vec<String>,
    },
    /// The schema attests but its singleton physical-ledger identity row is
    /// missing or malformed. Opening refuses rather than silently replacing
    /// the authority identity of an existing database.
    InstanceIdentityCorrupt {
        /// Bounded structural diagnosis.
        detail: String,
    },
    /// A new physical-ledger identity was required, but this build could not
    /// obtain operating-system entropy through its safe std-only backend.
    InstanceIdentityUnavailable {
        /// Entropy-source diagnosis.
        detail: String,
    },
    /// Identical bytes offered with a DIFFERENT envelope (kind or
    /// metadata) than the stored artifact carries (bead gp3.19).
    /// Byte deduplication stays content-addressed, but an envelope
    /// disagreement refuses instead of silently retaining whichever
    /// arrived first — provenance must not depend on insertion order.
    ArtifactEnvelopeConflict {
        /// Hex content hash of the artifact.
        hash_hex: String,
        /// Which envelope field conflicts: `kind` or `meta`.
        field: &'static str,
        /// The envelope value already stored (`<none>` for NULL meta).
        stored: String,
        /// The envelope value this call offered.
        offered: String,
    },
}

impl LedgerError {
    /// Stable machine-readable error code.
    #[must_use]
    pub fn code(&self) -> &'static str {
        match self {
            LedgerError::Open { .. } => "LedgerOpen",
            LedgerError::FutureSchema { .. } => "LedgerFutureSchema",
            LedgerError::Sql { .. } => "LedgerSql",
            LedgerError::Busy { .. } => "LedgerBusy",
            LedgerError::MissingExplicit { .. } => "LedgerMissingExplicit",
            LedgerError::Invalid { .. } => "LedgerInvalid",
            LedgerError::UnknownStateSemantics { .. } => "LedgerUnknownStateSemantics",
            LedgerError::Corrupt { .. } => "LedgerCorruption",
            LedgerError::ArtifactReadLimit { .. } => "LedgerArtifactReadLimit",
            LedgerError::OpCorrupt { .. } => "LedgerOpCorruption",
            LedgerError::TuneCorrupt { .. } => "LedgerTuneCorruption",
            LedgerError::TuneReadLimit { .. } => "LedgerTuneReadLimit",
            LedgerError::NotFound { .. } => "LedgerNotFound",
            LedgerError::DoubleFinish { .. } => "LedgerDoubleFinish",
            LedgerError::OpLineageSealed { .. } => "LedgerOpLineageSealed",
            LedgerError::WriterInTransaction => "LedgerWriterInTransaction",
            LedgerError::SchemaMismatch { .. } => "LedgerSchemaMismatch",
            LedgerError::InstanceIdentityCorrupt { .. } => "LedgerInstanceIdentityCorrupt",
            LedgerError::InstanceIdentityUnavailable { .. } => "LedgerInstanceIdentityUnavailable",
            LedgerError::ArtifactEnvelopeConflict { .. } => "LedgerArtifactEnvelopeConflict",
        }
    }
}

impl std::fmt::Display for LedgerError {
    #[allow(clippy::too_many_lines)] // One exhaustive agent-facing error vocabulary must stay visibly aligned with code().
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LedgerError::Open { path, detail } => {
                write!(f, "LedgerOpen: cannot open ledger at {path}: {detail}")
            }
            LedgerError::FutureSchema { found, supported } => write!(
                f,
                "LedgerFutureSchema: file has schema v{found} but this build supports up to \
                 v{supported} — upgrade fs-ledger; refusing to touch a newer file"
            ),
            LedgerError::Sql { context, detail } => {
                write!(f, "LedgerSql: {context}: {detail}")
            }
            LedgerError::Busy { context, detail } => write!(
                f,
                "LedgerBusy: {context}: {detail} — retryable contention; retry with backoff"
            ),
            LedgerError::MissingExplicit { field, problem } => write!(
                f,
                "LedgerMissingExplicit: Five Explicits field '{field}' rejected: {problem} \
                 (units travel inside the typed IR; seeds/versions/budget/capability are \
                 mandatory columns — Decalogue P4)"
            ),
            LedgerError::Invalid { field, problem } => {
                write!(f, "LedgerInvalid: field '{field}' rejected: {problem}")
            }
            LedgerError::UnknownStateSemantics {
                state_slot_hex,
                stored_law,
                stored_law_version,
                stored_state_schema_version,
                expected_law,
                expected_law_version,
                expected_state_schema_version,
                stored_parameters_hash,
                expected_parameters_hash,
                stored_contract_and_code_hash,
                expected_contract_and_code_hash,
                differences,
            } => write!(
                f,
                "LedgerUnknownStateSemantics: state slot {state_slot_hex} stores \
                 law={stored_law:?} law_version={stored_law_version} \
                 state_schema_version={stored_state_schema_version} \
                 parameters={stored_parameters_hash} implementation={stored_contract_and_code_hash}, \
                 but replay supplied law={expected_law:?} law_version={expected_law_version} \
                 state_schema_version={expected_state_schema_version} \
                 parameters={expected_parameters_hash} implementation={expected_contract_and_code_hash}; \
                 mismatched fields: {} — state bytes were withheld",
                differences.join(", ")
            ),
            LedgerError::Corrupt { hash_hex, detail } => write!(
                f,
                "LedgerCorruption: artifact {hash_hex}: {detail} — refuse to trust or replay \
                 from a tampered ledger"
            ),
            LedgerError::ArtifactReadLimit {
                hash_hex,
                limit,
                observed,
            } => write!(
                f,
                "LedgerArtifactReadLimit: stored metadata for artifact {hash_hex} declares \
                 {observed} bytes, exceeding the caller's {limit}-byte materialization budget"
            ),
            LedgerError::OpCorrupt { op, detail } => write!(
                f,
                "LedgerOpCorruption: op {op}: {detail} — refuse to materialize an operation \
                 row that bypassed the canonical bounded write path"
            ),
            LedgerError::TuneCorrupt { kernel, detail } => write!(
                f,
                "LedgerTuneCorruption: tune history for kernel {:?}: {detail} — refuse to \
                 materialize or interpret a row that bypassed the canonical tune API",
                kernel.get(..kernel.len().min(96)).unwrap_or(kernel)
            ),
            LedgerError::TuneReadLimit {
                kernel,
                resource,
                limit,
                observed_at_least,
            } => write!(
                f,
                "LedgerTuneReadLimit: tune history for kernel {:?} has {resource}=\
                 {observed_at_least}, exceeding limit {limit}; narrow or compact the history",
                kernel.get(..kernel.len().min(96)).unwrap_or(kernel)
            ),
            LedgerError::NotFound { what } => write!(f, "LedgerNotFound: {what}"),
            LedgerError::DoubleFinish { op } => write!(
                f,
                "LedgerDoubleFinish: op {op} already has an outcome; ops are event-sourced \
                 facts and cannot be finished twice"
            ),
            LedgerError::OpLineageSealed { op } => write!(
                f,
                "LedgerOpLineageSealed: op {op} already has an outcome; finish_op seals its \
                 input/output lineage, so record every edge before publishing the outcome"
            ),
            LedgerError::WriterInTransaction => write!(
                f,
                "LedgerWriterInTransaction: ArtifactWriter owns its own transaction; commit \
                 or roll back the explicit transaction first"
            ),
            LedgerError::SchemaMismatch {
                claimed_version,
                violations,
            } => write!(
                f,
                "LedgerSchemaMismatch: the file claims schema v{claimed_version} but its \
                 objects do not attest ({} violation(s)): {} — refusing to advance \
                 user_version or trust ledger authority through an alien or partially \
                 initialized schema; migrate the data out manually or delete the file if \
                 it is disposable",
                violations.len(),
                violations.join("; ")
            ),
            LedgerError::InstanceIdentityCorrupt { detail } => write!(
                f,
                "LedgerInstanceIdentityCorrupt: {detail}; refusing to replace or guess the \
                 physical ledger identity"
            ),
            LedgerError::InstanceIdentityUnavailable { detail } => write!(
                f,
                "LedgerInstanceIdentityUnavailable: {detail}; refusing to mint a predictable \
                 physical ledger identity"
            ),
            LedgerError::ArtifactEnvelopeConflict {
                hash_hex,
                field,
                stored,
                offered,
            } => write!(
                f,
                "LedgerArtifactEnvelopeConflict: artifact {hash_hex} already stores \
                 {field}={stored} but this call offered {field}={offered}; identical \
                 bytes dedupe only under an AGREEING envelope — match the stored \
                 envelope (or offer no metadata to accept it) instead of relying on \
                 insertion order",
            ),
        }
    }
}

impl std::error::Error for LedgerError {}

fn is_retryable(e: &FrankenError) -> bool {
    matches!(
        e,
        FrankenError::Busy
            | FrankenError::BusyRecovery
            | FrankenError::BusySnapshot { .. }
            | FrankenError::DatabaseLocked { .. }
            | FrankenError::WriteConflict { .. }
            | FrankenError::SerializationFailure { .. }
    )
}

fn is_duplicate_key(e: &FrankenError) -> bool {
    matches!(
        e,
        FrankenError::PrimaryKeyViolation | FrankenError::UniqueViolation { .. }
    )
}

fn sql_err(context: &str, e: &FrankenError) -> LedgerError {
    if is_retryable(e) {
        LedgerError::Busy {
            context: context.to_string(),
            detail: e.to_string(),
        }
    } else {
        LedgerError::Sql {
            context: context.to_string(),
            detail: e.to_string(),
        }
    }
}

// ---------------------------------------------------------------------------
// Public row/receipt types
// ---------------------------------------------------------------------------

/// Receipt returned by artifact writes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PutReceipt {
    /// Content identity.
    pub hash: ContentHash,
    /// Total byte length.
    pub len: u64,
    /// True if identical bytes were already stored (no new row).
    pub deduped: bool,
    /// True if stored as chunk rows rather than inline bytes.
    pub chunked: bool,
}

/// Metadata for one stored artifact.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ArtifactInfo {
    /// Content identity.
    pub hash: ContentHash,
    /// Caller-declared kind (e.g. "field", "mesh", "study-ir").
    pub kind: String,
    /// Total byte length.
    pub len: u64,
    /// Number of chunk rows (0 = stored inline).
    pub chunk_count: u64,
    /// Optional JSON metadata.
    pub meta: Option<String>,
    /// Wall-clock nanoseconds at first insertion (provenance envelope; not
    /// part of the content identity).
    pub created_at: i64,
}

/// Outcome of a finished op.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OpOutcome {
    /// Completed successfully.
    Ok,
    /// Failed with a structured diagnostic.
    Error,
    /// Cancelled (request → drain → finalize; P7).
    Cancelled,
}

impl OpOutcome {
    fn as_str(self) -> &'static str {
        match self {
            OpOutcome::Ok => "ok",
            OpOutcome::Error => "error",
            OpOutcome::Cancelled => "cancelled",
        }
    }
}

/// Lineage edge direction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeRole {
    /// The op consumed the artifact.
    In,
    /// The op produced the artifact.
    Out,
}

impl EdgeRole {
    fn as_str(self) -> &'static str {
        match self {
            EdgeRole::In => "in",
            EdgeRole::Out => "out",
        }
    }

    fn parse(value: &str) -> Option<Self> {
        match value {
            "in" => Some(Self::In),
            "out" => Some(Self::Out),
            _ => None,
        }
    }
}

/// Fixed-size execution context for one operation.
///
/// This intentionally excludes the variable-size IR and Five Explicits so a
/// verifier can check branch/mode provenance without materializing an
/// [`OpRow`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OpExecutionContext {
    /// Branch containing the operation.
    pub branch: i64,
    /// Recorded deterministic or fast execution mode.
    pub exec_mode: ExecMode,
}

/// One role-qualified artifact edge returned by a bounded lineage query.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OpArtifactEdge {
    /// Whether the operation consumed or produced the artifact.
    pub role: EdgeRole,
    /// Schema-v1 compatibility hash of the linked artifact. Use
    /// [`Ledger::edge_content_identity`] at typed-identity trust boundaries.
    pub artifact: ContentHash,
}

/// Capped output-producer lookup for one artifact.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BoundedProducerOps {
    /// Producer operation ids, ordered by id and limited to the caller cap.
    pub op_ids: Vec<i64>,
    /// `true` when at least one additional producer exists beyond `op_ids`.
    pub truncated: bool,
}

/// Capped role-qualified artifact edges for one operation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BoundedOpArtifactEdges {
    /// Edges ordered by role and artifact identity, limited to the caller cap.
    pub edges: Vec<OpArtifactEdge>,
    /// `true` when at least one additional edge exists beyond `edges`.
    pub truncated: bool,
}

/// The frozen Five Explicits of an op (P4). Units are the fifth explicit and
/// travel inside the typed IR itself (fs-qty dimensions), so they have no
/// separate column; the other four are mandatory here.
#[derive(Debug, Clone, Copy)]
pub struct FiveExplicits<'a> {
    /// RNG seed bytes (non-empty).
    pub seed: &'a [u8],
    /// Constellation/crate versions, JSON (e.g. the lock hash).
    pub versions: &'a str,
    /// Budget grant, JSON (accuracy/time/memory).
    pub budget: &'a str,
    /// Capability grant, JSON (ops/cores/mem/wall).
    pub capability: &'a str,
}

/// One recorded op row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OpRow {
    /// Rowid.
    pub id: i64,
    /// Session identity, if any.
    pub session: Option<Vec<u8>>,
    /// Frozen IR (JSON).
    pub ir: String,
    /// Frozen seed bytes.
    pub seed: Vec<u8>,
    /// Frozen versions (JSON).
    pub versions: String,
    /// Frozen budget (JSON).
    pub budget: String,
    /// Frozen capability (JSON).
    pub capability: String,
    /// Start wall-clock ns.
    pub t_start: i64,
    /// End wall-clock ns (None while in flight).
    pub t_end: Option<i64>,
    /// Outcome text (None while in flight).
    pub outcome: Option<String>,
    /// Structured diagnostic (JSON), if any.
    pub diag: Option<String>,
}

/// One event-stream row to append.
#[derive(Debug, Clone, Copy)]
pub struct EventRow<'a> {
    /// Session identity, if any.
    pub session: Option<&'a [u8]>,
    /// Logical or wall time (caller-controlled so deterministic replays can
    /// use logical time).
    pub t: i64,
    /// Event kind (fs-obs kind names recommended).
    pub kind: &'a str,
    /// JSON payload, if any.
    pub payload: Option<&'a str>,
}

/// One autotuner cache row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TuneRow {
    /// Kernel identity.
    pub kernel: String,
    /// Shape class.
    pub shape_class: String,
    /// Machine fingerprint bytes (fs-substrate probe hash).
    pub machine: Vec<u8>,
    /// Chosen parameters (JSON).
    pub params: String,
    /// Measured results (JSON).
    pub measured: String,
}

/// Rev S extension tables (sparse in v0; uniform `(name, body JSON)` shape).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExtensionTable {
    /// Requirement records.
    Requirements,
    /// Model cards for physical/surrogate models.
    ModelCards,
    /// Evidence records (Gauntlet artifacts, certificates).
    Evidence,
    /// Scenario definitions.
    Scenarios,
    /// Constraint definitions.
    Constraints,
    /// Hardware capability probe snapshots.
    CapabilityProbes,
    /// External import receipts.
    Imports,
    /// Unsafe-capsule registry mirror.
    UnsafeCapsules,
    /// Speculation telemetry: solve-node records carrying
    /// `(proposer_id, accepted, bound, iterations_saved)` (v3,
    /// bead lmp4.3).
    Speculation,
}

impl ExtensionTable {
    /// The underlying table name.
    #[must_use]
    pub fn table_name(self) -> &'static str {
        match self {
            ExtensionTable::Requirements => "requirements",
            ExtensionTable::Speculation => "speculation",
            ExtensionTable::ModelCards => "model_cards",
            ExtensionTable::Evidence => "evidence",
            ExtensionTable::Scenarios => "scenarios",
            ExtensionTable::Constraints => "constraints",
            ExtensionTable::CapabilityProbes => "capability_probes",
            ExtensionTable::Imports => "imports",
            ExtensionTable::UnsafeCapsules => "unsafe_capsules",
        }
    }
}

/// Referential/shape hygiene report. All-zero means the storage envelopes are
/// clean; semantic receipt replay and artifact re-hashing remain separate,
/// explicit verification operations.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct LintReport {
    /// Edges whose op id does not exist.
    pub orphan_edge_ops: u64,
    /// Edges whose artifact hash does not exist.
    pub orphan_edge_artifacts: u64,
    /// Metrics whose op id does not exist.
    pub orphan_metric_ops: u64,
    /// Artifacts violating the inline-XOR-chunked or bounded-row storage
    /// invariants.
    pub malformed_artifacts: u64,
    /// Chunked artifacts whose row count or dense zero-based sequence differs
    /// from metadata.
    pub chunk_count_mismatches: u64,
    /// Artifacts whose stored byte length differs from `len`.
    pub len_mismatches: u64,
    /// Chunk rows without a parent artifact row (e.g. abandoned staging).
    pub orphan_chunks: u64,
    /// Tune rows violating bounded storage types, lengths, or JSON validity.
    pub malformed_tune_rows: u64,
    /// Operation rows violating bounded storage types, lengths, JSON validity,
    /// or the start/finish outcome envelope.
    pub malformed_ops: u64,
    /// Ops with exactly one of (t_end, outcome) set.
    pub half_finished_ops: u64,
    /// Ops whose branch id does not exist (v2).
    pub orphan_op_branches: u64,
    /// Branches whose parent id does not exist (v2).
    pub orphan_branch_parents: u64,
    /// Quantity crosswalk rows whose historical source artifact is absent.
    pub orphan_qty_crosswalk_sources: u64,
    /// Quantity crosswalk rows whose canonical target artifact is absent.
    pub orphan_qty_crosswalk_targets: u64,
    /// Quantity crosswalk rows violating their fixed storage/version/rule
    /// envelope.
    pub malformed_qty_crosswalks: u64,
    /// Semantic state-checkpoint rows whose runtime-state artifact is absent.
    pub orphan_state_checkpoint_artifacts: u64,
    /// Semantic state-checkpoint rows violating their bounded SQL storage
    /// envelope. This does not re-hash receipts or runtime-state artifacts.
    pub malformed_state_checkpoints: u64,
}

impl LintReport {
    /// True when every counter is zero.
    #[must_use]
    pub fn is_clean(&self) -> bool {
        *self == LintReport::default()
    }
}

/// Result of a full integrity re-hash.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct IntegrityReport {
    /// Artifacts checked.
    pub checked: u64,
    /// Hex hashes of artifacts whose bytes no longer match their identity.
    pub corrupted: Vec<String>,
}

impl IntegrityReport {
    /// True when nothing is corrupted.
    #[must_use]
    pub fn is_clean(&self) -> bool {
        self.corrupted.is_empty()
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Wall-clock nanoseconds since the Unix epoch (provenance envelope only;
/// never part of content identity — P2).
#[must_use]
pub fn now_wall_ns() -> i64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| i64::try_from(d.as_nanos()).unwrap_or(i64::MAX))
}

/// Total i64 from a length (lengths are capped far below `i64::MAX` by the
/// engine's 1 GB value limit; saturate rather than wrap on absurd inputs).
fn int_from_usize(n: usize) -> i64 {
    i64::try_from(n).unwrap_or(i64::MAX)
}

fn int_from_u64(n: u64) -> i64 {
    i64::try_from(n).unwrap_or(i64::MAX)
}

fn normalize_schema_sql(sql: &str) -> String {
    sql.split_whitespace().collect::<Vec<_>>().join(" ")
}

const LEDGER_IDENTITY_GUARDS: [(&str, &str); 3] = [
    ("trg_ledger_identity_immutable_update", schema::V5[0]),
    ("trg_ledger_identity_immutable_delete", schema::V5[1]),
    ("trg_ledger_identity_immutable_reinsert", schema::V5[2]),
];

/// Every user object in `sqlite_master` as `(type, name) -> normalized
/// SQL` (whitespace-collapsed; internal `sqlite_*` entries and DDL-less
/// auto-indexes excluded). The stored SQL text carries STRICT, CHECKs,
/// and FOREIGN KEY clauses verbatim, so equality attests them all.
fn schema_objects(
    conn: &Connection,
) -> Result<std::collections::BTreeMap<(String, String), String>, LedgerError> {
    let rows = conn
        .query("SELECT type, name, COALESCE(sql, '') FROM sqlite_master ORDER BY type, name")
        .map_err(|e| sql_err("attest: read sqlite_master", &e))?;
    let mut objects = std::collections::BTreeMap::new();
    for row in &rows {
        let kind = row_text(row, 0, "attest: sqlite_master type")?;
        let name = row_text(row, 1, "attest: sqlite_master name")?;
        // SQL LIKE treats `_` as a wildcard, so filtering with
        // `NOT LIKE 'sqlite_%'` would also hide legal user objects such as
        // `sqlitex_hidden`. Only the literal reserved prefix is internal.
        if name.starts_with("sqlite_") {
            continue;
        }
        let sql = row_text(row, 2, "attest: sqlite_master sql")?;
        let normalized = normalize_schema_sql(&sql);
        objects.insert((kind, name), normalized);
    }
    Ok(objects)
}

/// One table's columns as stable signature strings
/// (`name:type:notnull:default:pk`) from `PRAGMA table_info`.
fn table_columns(conn: &Connection, table: &str) -> Result<Vec<String>, LedgerError> {
    let rows = conn
        .query(&format!("PRAGMA table_info({table})"))
        .map_err(|e| sql_err("attest: table_info", &e))?;
    let mut columns = Vec::with_capacity(rows.len());
    for row in &rows {
        let name = row_text(row, 1, "attest: column name")?;
        let declared = row_text(row, 2, "attest: column type")?;
        let not_null = row_i64(row, 3, "attest: column not-null")? != 0;
        let default_sql = match row.get(4) {
            Some(SqliteValue::Text(value)) => value.as_str().to_string(),
            Some(SqliteValue::Null) | None => "<none>".to_string(),
            Some(other) => format!("{other:?}"),
        };
        let pk = row_i64(row, 5, "attest: column pk")? != 0;
        columns.push(format!("{name}:{declared}:{not_null}:{default_sql}:{pk}"));
    }
    columns.sort();
    Ok(columns)
}

/// A fresh in-memory database with the first `steps` shipped migration
/// batches applied — the attestation reference.
fn reference_connection(steps: usize) -> Result<Connection, LedgerError> {
    let reference = Connection::open(":memory:").map_err(|e| LedgerError::Sql {
        context: "attest: open reference".to_string(),
        detail: e.to_string(),
    })?;
    for batch in schema::MIGRATIONS.iter().take(steps) {
        for ddl in *batch {
            reference.execute(ddl).map_err(|e| LedgerError::Sql {
                context: "attest: build reference".to_string(),
                detail: format!("{e} while executing: {}", ddl.get(..60).unwrap_or(ddl)),
            })?;
        }
    }
    Ok(reference)
}

fn row_text(row: &fsqlite::Row, idx: usize, context: &str) -> Result<String, LedgerError> {
    match row.get(idx) {
        Some(SqliteValue::Text(value)) => Ok(value.as_str().to_string()),
        other => Err(LedgerError::Sql {
            context: context.to_string(),
            detail: format!("expected TEXT at column {idx}, got {other:?}"),
        }),
    }
}

fn text_param(s: &str) -> SqliteValue {
    SqliteValue::Text(s.into())
}

fn blob_param(b: &[u8]) -> SqliteValue {
    SqliteValue::Blob(b.to_vec().into())
}

fn opt_text_param(s: Option<&str>) -> SqliteValue {
    s.map_or(SqliteValue::Null, text_param)
}

fn opt_blob_param(b: Option<&[u8]>) -> SqliteValue {
    b.map_or(SqliteValue::Null, blob_param)
}

fn row_i64(row: &fsqlite::Row, idx: usize, context: &str) -> Result<i64, LedgerError> {
    match row.get(idx) {
        Some(SqliteValue::Integer(v)) => Ok(*v),
        other => Err(LedgerError::Sql {
            context: context.to_string(),
            detail: format!("column {idx}: expected INTEGER, got {other:?}"),
        }),
    }
}

fn nonnegative_u64(value: i64, context: &str) -> Result<u64, LedgerError> {
    u64::try_from(value).map_err(|_| LedgerError::Sql {
        context: context.to_string(),
        detail: format!("expected non-negative INTEGER, got {value}"),
    })
}

fn row_u64(row: &fsqlite::Row, idx: usize, context: &str) -> Result<u64, LedgerError> {
    nonnegative_u64(row_i64(row, idx, context)?, context)
}

fn validate_tune_identity(field: &str, value: &str, max_bytes: usize) -> Result<(), LedgerError> {
    if value.is_empty() {
        return Err(LedgerError::Invalid {
            field: field.to_string(),
            problem: "empty; tune identities must be non-empty visible ASCII".to_string(),
        });
    }
    if value.len() > max_bytes {
        return Err(LedgerError::Invalid {
            field: field.to_string(),
            problem: format!(
                "{} bytes exceeds the {max_bytes}-byte tune identity limit",
                value.len()
            ),
        });
    }
    if let Some((offset, byte)) = value
        .bytes()
        .enumerate()
        .find(|(_, byte)| !(b'!'..=b'~').contains(byte))
    {
        return Err(LedgerError::Invalid {
            field: field.to_string(),
            problem: format!(
                "byte 0x{byte:02x} at offset {offset} is not visible ASCII; use bytes 0x21..=0x7e"
            ),
        });
    }
    Ok(())
}

fn validate_tune_machine(machine: &[u8]) -> Result<(), LedgerError> {
    if machine.is_empty() {
        return Err(LedgerError::Invalid {
            field: "machine".to_string(),
            problem: "empty; supply an exact machine fingerprint blob".to_string(),
        });
    }
    if machine.len() > MAX_TUNE_MACHINE_BYTES {
        return Err(LedgerError::Invalid {
            field: "machine".to_string(),
            problem: format!(
                "{} bytes exceeds the {MAX_TUNE_MACHINE_BYTES}-byte machine fingerprint limit",
                machine.len()
            ),
        });
    }
    Ok(())
}

fn tune_corrupt(kernel: &str, detail: impl Into<String>) -> LedgerError {
    LedgerError::TuneCorrupt {
        kernel: kernel.to_string(),
        detail: detail.into(),
    }
}

fn op_storage_predicate() -> String {
    format!(
        "(session IS NULL OR (typeof(session) = 'blob' AND length(session) <= {MAX_OP_SESSION_BYTES})) AND \
         CASE WHEN typeof(ir) = 'text' THEN \
             CASE WHEN length(CAST(ir AS BLOB)) BETWEEN 1 AND {MAX_OP_IR_BYTES} \
                  THEN json_valid(ir) ELSE 0 END ELSE 0 END = 1 AND \
         typeof(seed) = 'blob' AND length(seed) BETWEEN 1 AND {MAX_OP_SEED_BYTES} AND \
         CASE WHEN typeof(versions) = 'text' THEN \
             CASE WHEN length(CAST(versions AS BLOB)) BETWEEN 1 AND {MAX_OP_VERSIONS_BYTES} \
                  THEN json_valid(versions) ELSE 0 END ELSE 0 END = 1 AND \
         CASE WHEN typeof(budget) = 'text' THEN \
             CASE WHEN length(CAST(budget AS BLOB)) BETWEEN 1 AND {MAX_OP_BUDGET_BYTES} \
                  THEN json_valid(budget) ELSE 0 END ELSE 0 END = 1 AND \
         CASE WHEN typeof(capability) = 'text' THEN \
             CASE WHEN length(CAST(capability AS BLOB)) BETWEEN 1 AND {MAX_OP_CAPABILITY_BYTES} \
                  THEN json_valid(capability) ELSE 0 END ELSE 0 END = 1 AND \
         typeof(t_start) = 'integer' AND \
         typeof(branch) = 'integer' AND \
         ((t_end IS NULL AND outcome IS NULL) OR \
          (typeof(t_end) = 'integer' AND typeof(outcome) = 'text' AND \
           outcome IN ('ok','error','cancelled'))) AND \
         typeof(exec_mode) = 'text' AND exec_mode IN ('deterministic','fast') AND \
         CASE WHEN diag IS NULL THEN 1 WHEN typeof(diag) = 'text' THEN \
             CASE WHEN length(CAST(diag AS BLOB)) BETWEEN 1 AND {MAX_OP_DIAG_BYTES} \
                  THEN json_valid(diag) ELSE 0 END ELSE 0 END = 1"
    )
}

fn bounded_lineage_query_limit(cap: usize) -> Result<usize, LedgerError> {
    if cap > MAX_LINEAGE_QUERY_ROWS {
        return Err(LedgerError::Invalid {
            field: "cap".to_string(),
            problem: format!(
                "lineage query cap {cap} exceeds the public maximum \
                 {MAX_LINEAGE_QUERY_ROWS}; narrow the verification query"
            ),
        });
    }
    cap.checked_add(1).ok_or_else(|| LedgerError::Invalid {
        field: "cap".to_string(),
        problem: "lineage query cap cannot be represented as cap + 1".to_string(),
    })
}

#[derive(Debug, Clone, Copy)]
struct TuneColumnSpec {
    type_idx: usize,
    len_idx: usize,
    field: &'static str,
    expected_type: &'static str,
    min_bytes: usize,
    max_bytes: usize,
}

fn tune_column_len(
    row: &fsqlite::Row,
    spec: TuneColumnSpec,
    kernel: &str,
) -> Result<usize, LedgerError> {
    let stored_type = row_text(row, spec.type_idx, "tune metadata type")?;
    if stored_type != spec.expected_type {
        return Err(tune_corrupt(
            kernel,
            format!(
                "{} has storage type {stored_type:?}, expected {}",
                spec.field, spec.expected_type
            ),
        ));
    }
    let raw_len = row_i64(row, spec.len_idx, "tune metadata length")?;
    let len = usize::try_from(raw_len).map_err(|_| {
        tune_corrupt(
            kernel,
            format!(
                "{} has negative or unrepresentable byte length {raw_len}",
                spec.field
            ),
        )
    })?;
    if len < spec.min_bytes || len > spec.max_bytes {
        return Err(tune_corrupt(
            kernel,
            format!(
                "{} byte length {len} is outside {}..={}",
                spec.field, spec.min_bytes, spec.max_bytes
            ),
        ));
    }
    Ok(len)
}

fn require_stored_json(
    row: &fsqlite::Row,
    valid_idx: usize,
    field: &str,
    kernel: &str,
) -> Result<(), LedgerError> {
    if row_i64(row, valid_idx, "tune metadata json_valid")? != 1 {
        return Err(tune_corrupt(kernel, format!("{field} is not valid JSON")));
    }
    Ok(())
}

fn stored_tune_identity(
    kernel: &str,
    field: &str,
    value: &str,
    max_bytes: usize,
) -> Result<(), LedgerError> {
    validate_tune_identity(field, value, max_bytes).map_err(|error| match error {
        LedgerError::Invalid { problem, .. } => tune_corrupt(kernel, format!("{field}: {problem}")),
        other => other,
    })
}

#[derive(Debug, Clone, Copy)]
struct TuneRowPreflight {
    materialized_bytes: usize,
}

fn tune_row_preflight(row: &fsqlite::Row, kernel: &str) -> Result<TuneRowPreflight, LedgerError> {
    let shape_len = tune_column_len(
        row,
        TuneColumnSpec {
            type_idx: 0,
            len_idx: 1,
            field: "shape_class",
            expected_type: "text",
            min_bytes: 1,
            max_bytes: MAX_TUNE_SHAPE_CLASS_BYTES,
        },
        kernel,
    )?;
    if row_i64(row, 2, "tune metadata canonical shape")? != 1 {
        return Err(tune_corrupt(
            kernel,
            "shape_class is not canonical visible ASCII",
        ));
    }
    let machine_len = tune_column_len(
        row,
        TuneColumnSpec {
            type_idx: 3,
            len_idx: 4,
            field: "machine",
            expected_type: "blob",
            min_bytes: 1,
            max_bytes: MAX_TUNE_MACHINE_BYTES,
        },
        kernel,
    )?;
    let params_len = tune_column_len(
        row,
        TuneColumnSpec {
            type_idx: 5,
            len_idx: 6,
            field: "params",
            expected_type: "text",
            min_bytes: 1,
            max_bytes: MAX_TUNE_PARAMS_BYTES,
        },
        kernel,
    )?;
    require_stored_json(row, 7, "params", kernel)?;
    let measured_len = tune_column_len(
        row,
        TuneColumnSpec {
            type_idx: 8,
            len_idx: 9,
            field: "measured",
            expected_type: "text",
            min_bytes: 1,
            max_bytes: MAX_TUNE_MEASURED_BYTES,
        },
        kernel,
    )?;
    require_stored_json(row, 10, "measured", kernel)?;
    let materialized_bytes = kernel
        .len()
        .checked_add(shape_len)
        .and_then(|sum| sum.checked_add(machine_len))
        .and_then(|sum| sum.checked_add(params_len))
        .and_then(|sum| sum.checked_add(measured_len))
        .ok_or_else(|| LedgerError::TuneReadLimit {
            kernel: kernel.to_string(),
            resource: "materialized_bytes",
            limit: MAX_TUNE_SCAN_BYTES,
            observed_at_least: usize::MAX,
        })?;
    Ok(TuneRowPreflight { materialized_bytes })
}

fn preflight_tune_scan(conn: &Connection, kernel: &str) -> Result<usize, LedgerError> {
    let metadata_sql = format!(
        "SELECT typeof(shape_class), length(CAST(shape_class AS BLOB)), \
                length(CAST(shape_class AS BLOB)) = length(shape_class) AND \
                    shape_class NOT GLOB '*[^!-~]*', \
                typeof(machine), length(machine), \
                typeof(params), length(CAST(params AS BLOB)), \
                    CASE WHEN typeof(params) = 'text' THEN \
                        CASE WHEN length(CAST(params AS BLOB)) BETWEEN 1 AND {MAX_TUNE_PARAMS_BYTES} \
                            THEN json_valid(params) ELSE 0 END \
                    ELSE 0 END, \
                typeof(measured), length(CAST(measured AS BLOB)), \
                    CASE WHEN typeof(measured) = 'text' THEN \
                        CASE WHEN length(CAST(measured AS BLOB)) BETWEEN 1 AND {MAX_TUNE_MEASURED_BYTES} \
                            THEN json_valid(measured) ELSE 0 END \
                    ELSE 0 END \
         FROM tune WHERE kernel = ?1 LIMIT {}",
        MAX_TUNE_ROWS_PER_KERNEL + 1
    );
    let metadata = conn
        .query_with_params(&metadata_sql, &[text_param(kernel)])
        .map_err(|e| sql_err("tune scan metadata", &e))?;
    if metadata.len() > MAX_TUNE_ROWS_PER_KERNEL {
        return Err(LedgerError::TuneReadLimit {
            kernel: kernel.to_string(),
            resource: "rows",
            limit: MAX_TUNE_ROWS_PER_KERNEL,
            observed_at_least: MAX_TUNE_ROWS_PER_KERNEL + 1,
        });
    }
    let mut aggregate_bytes = 0_usize;
    for row in &metadata {
        let preflight = tune_row_preflight(row, kernel)?;
        aggregate_bytes = aggregate_bytes
            .checked_add(preflight.materialized_bytes)
            .ok_or_else(|| LedgerError::TuneReadLimit {
                kernel: kernel.to_string(),
                resource: "materialized_bytes",
                limit: MAX_TUNE_SCAN_BYTES,
                observed_at_least: usize::MAX,
            })?;
        if aggregate_bytes > MAX_TUNE_SCAN_BYTES {
            return Err(LedgerError::TuneReadLimit {
                kernel: kernel.to_string(),
                resource: "materialized_bytes",
                limit: MAX_TUNE_SCAN_BYTES,
                observed_at_least: aggregate_bytes,
            });
        }
    }
    Ok(metadata.len())
}

fn guarded_tune_scan_sql(kernel_len: usize) -> String {
    format!(
        "SELECT shape_class, machine, params, measured FROM tune \
         WHERE kernel = ?1 AND \
               (SELECT COUNT(*) FROM tune AS tune_count \
                WHERE tune_count.kernel = ?1) <= {MAX_TUNE_ROWS_PER_KERNEL} AND \
               (SELECT COALESCE(SUM(\
                    length(CAST(tune_budget.shape_class AS BLOB)) + \
                    length(tune_budget.machine) + \
                    length(CAST(tune_budget.params AS BLOB)) + \
                    length(CAST(tune_budget.measured AS BLOB)) + {kernel_len}\
                ), 0) FROM tune AS tune_budget \
                WHERE tune_budget.kernel = ?1) <= {MAX_TUNE_SCAN_BYTES} AND \
               NOT EXISTS (SELECT 1 FROM tune AS tune_guard \
                   WHERE tune_guard.kernel = ?1 AND (\
                       typeof(tune_guard.shape_class) != 'text' OR \
                       length(CAST(tune_guard.shape_class AS BLOB)) NOT BETWEEN 1 AND {MAX_TUNE_SHAPE_CLASS_BYTES} OR \
                       length(CAST(tune_guard.shape_class AS BLOB)) != length(tune_guard.shape_class) OR \
                       tune_guard.shape_class GLOB '*[^!-~]*' OR \
                       typeof(tune_guard.machine) != 'blob' OR \
                       length(tune_guard.machine) NOT BETWEEN 1 AND {MAX_TUNE_MACHINE_BYTES} OR \
                       typeof(tune_guard.params) != 'text' OR \
                       length(CAST(tune_guard.params AS BLOB)) NOT BETWEEN 1 AND {MAX_TUNE_PARAMS_BYTES} OR \
                       CASE WHEN typeof(tune_guard.params) = 'text' THEN \
                           CASE WHEN length(CAST(tune_guard.params AS BLOB)) BETWEEN 1 AND {MAX_TUNE_PARAMS_BYTES} \
                               THEN json_valid(tune_guard.params) ELSE 0 END \
                           ELSE 0 END != 1 OR \
                       typeof(tune_guard.measured) != 'text' OR \
                       length(CAST(tune_guard.measured AS BLOB)) NOT BETWEEN 1 AND {MAX_TUNE_MEASURED_BYTES} OR \
                       CASE WHEN typeof(tune_guard.measured) = 'text' THEN \
                           CASE WHEN length(CAST(tune_guard.measured AS BLOB)) BETWEEN 1 AND {MAX_TUNE_MEASURED_BYTES} \
                               THEN json_valid(tune_guard.measured) ELSE 0 END \
                           ELSE 0 END != 1)) \
         ORDER BY shape_class, machine LIMIT {}",
        MAX_TUNE_ROWS_PER_KERNEL + 1
    )
}

fn materialize_tune_rows(kernel: &str, rows: &[fsqlite::Row]) -> Result<Vec<TuneRow>, LedgerError> {
    let mut out = Vec::with_capacity(rows.len());
    let mut materialized_bytes = 0_usize;
    for row in rows {
        let shape_class = row_text(row, 0, "tune guarded scan shape_class")?;
        stored_tune_identity(
            kernel,
            "shape_class",
            &shape_class,
            MAX_TUNE_SHAPE_CLASS_BYTES,
        )?;
        let machine = match row.get(1) {
            Some(SqliteValue::Blob(bytes)) => bytes.to_vec(),
            other => {
                return Err(tune_corrupt(
                    kernel,
                    format!("machine: expected BLOB, got {other:?}"),
                ));
            }
        };
        validate_tune_machine(&machine).map_err(|error| match error {
            LedgerError::Invalid { problem, .. } => {
                tune_corrupt(kernel, format!("machine: {problem}"))
            }
            other => other,
        })?;
        let params = row_text(row, 2, "tune guarded scan params")?;
        let measured = row_text(row, 3, "tune guarded scan measured")?;
        materialized_bytes = materialized_bytes
            .checked_add(kernel.len())
            .and_then(|sum| sum.checked_add(shape_class.len()))
            .and_then(|sum| sum.checked_add(machine.len()))
            .and_then(|sum| sum.checked_add(params.len()))
            .and_then(|sum| sum.checked_add(measured.len()))
            .ok_or_else(|| LedgerError::TuneReadLimit {
                kernel: kernel.to_string(),
                resource: "materialized_bytes",
                limit: MAX_TUNE_SCAN_BYTES,
                observed_at_least: usize::MAX,
            })?;
        out.push(TuneRow {
            kernel: kernel.to_string(),
            shape_class,
            machine,
            params,
            measured,
        });
    }
    if materialized_bytes > MAX_TUNE_SCAN_BYTES {
        return Err(LedgerError::TuneReadLimit {
            kernel: kernel.to_string(),
            resource: "materialized_bytes",
            limit: MAX_TUNE_SCAN_BYTES,
            observed_at_least: materialized_bytes,
        });
    }
    Ok(out)
}

fn guarded_tune_scan(
    conn: &Connection,
    kernel: &str,
    expected_rows: usize,
) -> Result<Vec<TuneRow>, LedgerError> {
    let rows = conn
        .query_with_params(&guarded_tune_scan_sql(kernel.len()), &[text_param(kernel)])
        .map_err(|e| sql_err("tune guarded scan", &e))?;
    if rows.len() != expected_rows {
        return Err(LedgerError::Busy {
            context: "tune guarded scan".to_string(),
            detail: format!(
                "history changed after bounded metadata preflight ({expected_rows} rows became {})",
                rows.len()
            ),
        });
    }
    materialize_tune_rows(kernel, &rows)
}

fn inline_artifact_bytes(row: &fsqlite::Row, expected_len: u64) -> Result<&[u8], String> {
    let Some(SqliteValue::Blob(bytes)) = row.get(0) else {
        return Err(format!("inline bytes: expected BLOB, got {:?}", row.get(0)));
    };
    let actual_len = u64::try_from(bytes.len())
        .map_err(|_| "inline byte length does not fit the ledger length domain".to_string())?;
    if actual_len != expected_len {
        return Err(format!(
            "inline length mismatch: recorded {expected_len}, found {actual_len}"
        ));
    }
    Ok(bytes.as_ref())
}

fn storage_i64(row: &fsqlite::Row, idx: usize, field: &str) -> Result<i64, String> {
    match row.get(idx) {
        Some(SqliteValue::Integer(value)) => Ok(*value),
        other => Err(format!(
            "storage preflight {field}: expected INTEGER at column {idx}, got {other:?}"
        )),
    }
}

fn storage_u64(row: &fsqlite::Row, idx: usize, field: &str) -> Result<u64, String> {
    let value = storage_i64(row, idx, field)?;
    u64::try_from(value)
        .map_err(|_| format!("storage preflight {field}: expected non-negative value, got {value}"))
}

#[derive(Debug, Clone, Copy)]
struct ArtifactChunkPreflight {
    count: u64,
    non_null_count: u64,
    min_seq: i64,
    max_seq: i64,
    total_len: u64,
    max_len: u64,
}

enum ChunkedArtifactInsert {
    Inserted,
    Deduped(ArtifactInfo),
}

impl ArtifactChunkPreflight {
    fn from_row(row: &fsqlite::Row) -> Result<Self, String> {
        Ok(Self {
            count: storage_u64(row, 0, "chunk count")?,
            non_null_count: storage_u64(row, 1, "non-null chunk count")?,
            min_seq: storage_i64(row, 2, "minimum chunk sequence")?,
            max_seq: storage_i64(row, 3, "maximum chunk sequence")?,
            total_len: storage_u64(row, 4, "chunk byte total")?,
            max_len: storage_u64(row, 5, "maximum chunk length")?,
        })
    }

    fn validate(self, info: &ArtifactInfo) -> Result<(), String> {
        let storage_chunk_len = STORAGE_CHUNK_LEN as u64;
        if self.non_null_count != self.count {
            return Err(format!(
                "chunk BLOB count mismatch: found {} rows but only {} non-NULL byte values",
                self.count, self.non_null_count
            ));
        }
        if self.max_len > storage_chunk_len {
            return Err(format!(
                "chunk exceeds the {STORAGE_CHUNK_LEN}-byte storage bound: maximum stored length is {} bytes",
                self.max_len
            ));
        }
        if self.count != info.chunk_count {
            return Err(format!(
                "chunk count mismatch: recorded {}, found {}",
                info.chunk_count, self.count
            ));
        }
        let expected_max = i64::try_from(info.chunk_count - 1)
            .map_err(|_| "recorded chunk count does not fit the sequence domain".to_string())?;
        if self.min_seq != 0 || self.max_seq != expected_max {
            return Err(format!(
                "chunk sequence range is not dense from zero: count {}, min {}, max {}",
                self.count, self.min_seq, self.max_seq
            ));
        }
        if self.total_len != info.len {
            return Err(format!(
                "chunk length mismatch: recorded {}, found {}",
                info.len, self.total_len
            ));
        }
        Ok(())
    }
}

struct ArtifactChunkValidator {
    expected_len: u64,
    expected_count: u64,
    streamed: u64,
    actual_count: u64,
}

impl ArtifactChunkValidator {
    const fn new(info: &ArtifactInfo) -> Self {
        Self {
            expected_len: info.len,
            expected_count: info.chunk_count,
            streamed: 0,
            actual_count: 0,
        }
    }

    fn accept<'row>(&mut self, row: &'row fsqlite::Row) -> Result<&'row [u8], String> {
        let actual_seq = match row.get(0) {
            Some(SqliteValue::Integer(seq)) => u64::try_from(*seq)
                .map_err(|_| format!("chunk sequence must be non-negative, got {seq}"))?,
            other => {
                return Err(format!("chunk sequence: expected INTEGER, got {other:?}"));
            }
        };
        if actual_seq != self.actual_count {
            return Err(format!(
                "chunk sequence mismatch: expected {}, found {actual_seq}",
                self.actual_count
            ));
        }

        let Some(SqliteValue::Blob(bytes)) = row.get(1) else {
            return Err(format!("chunk bytes: expected BLOB, got {:?}", row.get(1)));
        };
        if bytes.len() > STORAGE_CHUNK_LEN {
            return Err(format!(
                "chunk {actual_seq} exceeds the {STORAGE_CHUNK_LEN}-byte storage bound: found {} \
                 bytes",
                bytes.len()
            ));
        }
        let chunk_len = u64::try_from(bytes.len()).map_err(|_| {
            format!("chunk {actual_seq} length does not fit the ledger length domain")
        })?;
        let next_count = self
            .actual_count
            .checked_add(1)
            .ok_or_else(|| "actual chunk count overflowed u64".to_string())?;
        if next_count > self.expected_count {
            return Err(format!(
                "chunk count exceeds recorded {} at sequence {actual_seq}",
                self.expected_count
            ));
        }
        let next_streamed = self
            .streamed
            .checked_add(chunk_len)
            .ok_or_else(|| format!("byte total overflowed u64 at chunk sequence {actual_seq}"))?;
        if next_streamed > self.expected_len {
            return Err(format!(
                "chunk bytes exceed recorded length {} at sequence {actual_seq}",
                self.expected_len
            ));
        }

        self.actual_count = next_count;
        self.streamed = next_streamed;
        Ok(bytes.as_ref())
    }

    fn finish(self) -> Result<u64, String> {
        if self.actual_count != self.expected_count {
            return Err(format!(
                "chunk count mismatch: recorded {}, found {}",
                self.expected_count, self.actual_count
            ));
        }
        if self.streamed != self.expected_len {
            return Err(format!(
                "chunk length mismatch: recorded {}, streamed {}",
                self.expected_len, self.streamed
            ));
        }
        Ok(self.streamed)
    }
}

// ---------------------------------------------------------------------------
// The Ledger
// ---------------------------------------------------------------------------

/// One handle on the Design Ledger: a single fsqlite connection plus the
/// schema/pragma contract. Not `Send`: open one per thread (the engine's
/// documented model); snapshot isolation happens below this API.
pub struct Ledger {
    conn: Connection,
    path: String,
    instance_id: LedgerInstanceId,
    /// Monotone count of read-side queries issued through the typed read
    /// APIs (bead vm3i): the measurable basis for verification query
    /// budgets. Diagnostic only — never part of any receipt or identity.
    read_queries: core::cell::Cell<u64>,
}

impl Ledger {
    /// Open (creating and migrating as needed) the ledger at `path`.
    ///
    /// Applies the pragma contract: WAL journal, `synchronous=FULL`
    /// (fsync-before-publish durability), `busy_timeout`, and enforced
    /// foreign keys. `":memory:"` is supported for tests.
    ///
    /// # Errors
    /// [`LedgerError::Open`] on engine failure; [`LedgerError::FutureSchema`]
    /// if the file was written by a newer fs-ledger.
    pub fn open(path: &str) -> Result<Ledger, LedgerError> {
        let conn = Connection::open(path).map_err(|e| LedgerError::Open {
            path: path.to_string(),
            detail: e.to_string(),
        })?;
        for pragma in [
            "PRAGMA journal_mode=WAL",
            "PRAGMA synchronous=FULL",
            "PRAGMA busy_timeout=5000",
            "PRAGMA foreign_keys=ON",
        ] {
            conn.query(pragma).map_err(|e| LedgerError::Open {
                path: path.to_string(),
                detail: format!("{pragma}: {e}"),
            })?;
        }
        let mut ledger = Ledger {
            conn,
            path: path.to_string(),
            // Never escapes `open`: migration/readback below replaces this
            // sentinel with the persisted, shape-checked identity.
            instance_id: LedgerInstanceId([0; 16]),
            read_queries: core::cell::Cell::new(0),
        };
        ledger.migrate()?;
        ledger.instance_id = ledger.read_current_instance_id()?;
        Ok(ledger)
    }

    /// Monotone count of typed read-API queries issued through this
    /// connection (bead vm3i): `tune_rows`/`tune_get`/`tune_content_identity`,
    /// `get_artifact`/`op`, bounded lineage reads, seal reads, `edge_exists`,
    /// and `checked_instance_id` each count once. The measurable basis for
    /// verification query budgets; diagnostic only, never part of a receipt.
    #[must_use]
    pub fn read_queries(&self) -> u64 {
        self.read_queries.get()
    }

    fn note_read_query(&self) {
        self.read_queries
            .set(self.read_queries.get().saturating_add(1));
    }

    /// The path this ledger was opened at.
    #[must_use]
    pub fn path(&self) -> &str {
        &self.path
    }

    /// Cached open-time identity of this physical ledger instance.
    ///
    /// This is the authority key for binding higher-level sinks. It is stable
    /// across aliases and reopenings of one file, distinct for a replacement
    /// file at the same path, and distinct for independent in-memory handles.
    /// Use [`Ledger::checked_instance_id`] at a trust boundary that must also
    /// detect out-of-band mutation after this handle opened.
    #[must_use]
    pub const fn instance_id(&self) -> LedgerInstanceId {
        self.instance_id
    }

    /// Re-read and validate the persisted identity, then compare it with this
    /// handle's cached open-time authority.
    ///
    /// This is the fail-closed accessor for persistence and audit boundaries.
    /// It cheaply re-attests the three exact v5 identity-guard definitions
    /// before reading the row. The comparison also detects a client that
    /// bypassed and then restored those triggers after changing the identity
    /// while this handle remained open.
    ///
    /// # Errors
    /// [`LedgerError::SchemaMismatch`] if an identity guard is missing or
    /// differs from the shipped definition; [`LedgerError::InstanceIdentityCorrupt`]
    /// if the row is missing, malformed, or differs from the cached open-time
    /// identity. Engine errors are returned as [`LedgerError::Sql`] or
    /// [`LedgerError::Busy`].
    pub fn checked_instance_id(&self) -> Result<LedgerInstanceId, LedgerError> {
        self.note_read_query();
        self.attest_identity_guards()?;
        let current = self.read_current_instance_id()?;
        if current != self.instance_id {
            return Err(LedgerError::InstanceIdentityCorrupt {
                detail: format!(
                    "open handle cached identity {} but the current row contains {}",
                    self.instance_id, current
                ),
            });
        }
        Ok(current)
    }

    /// The schema version recorded in the file.
    ///
    /// # Errors
    /// [`LedgerError::Sql`] on engine failure.
    pub fn schema_version(&self) -> Result<i64, LedgerError> {
        let row = self
            .conn
            .query_row("PRAGMA user_version")
            .map_err(|e| sql_err("read user_version", &e))?;
        row_i64(&row, 0, "read user_version")
    }

    fn migrate(&self) -> Result<(), LedgerError> {
        let found = self.schema_version()?;
        self.migrate_from_observed_version(found)
    }

    fn migrate_from_observed_version(&self, found: i64) -> Result<(), LedgerError> {
        if found > SCHEMA_VERSION {
            return Err(LedgerError::FutureSchema {
                found,
                supported: SCHEMA_VERSION,
            });
        }
        // ATTESTATION BEFORE ADVANCEMENT (bead gp3.18): `CREATE TABLE IF
        // NOT EXISTS` silently tolerates pre-existing objects, so without
        // this check an alien or partially initialized file would be
        // stamped as the current schema. A v0 file must be EMPTY; a file
        // claiming v>0 must attest object-for-object against a reference
        // built from the shipped DDL for that version.
        if found == 0 {
            // ATOMIC INITIALIZATION: the whole ladder plus the final
            // version stamp in ONE transaction. The emptiness attestation is
            // inside that same transaction so another connection cannot add
            // an object between the check and the first CREATE statement.
            self.conn
                .begin_transaction()
                .map_err(|e| sql_err("init: begin", &e))?;
            // The preflight marker may be stale: another opener can finish
            // legitimate first initialization after our v0 read but before
            // this transaction begins. Refresh under the transaction before
            // treating any objects as alien v0 state, then leave the stale
            // snapshot and run the ordinary attestation/migration path.
            let found_after_begin = match self.schema_version() {
                Ok(found_after_begin) => found_after_begin,
                Err(error) => {
                    let _ = self.conn.rollback_transaction();
                    return Err(error);
                }
            };
            if found_after_begin != 0 {
                self.conn
                    .rollback_transaction()
                    .map_err(|error| sql_err("init: roll back stale v0 observation", &error))?;
                return self.migrate();
            }
            let init = (|| {
                let objects = schema_objects(&self.conn)?;
                if !objects.is_empty() {
                    let violations = objects
                        .keys()
                        .map(|(kind, name)| {
                            format!("pre-existing {kind} `{name}` in an unversioned file")
                        })
                        .collect();
                    return Err(LedgerError::SchemaMismatch {
                        claimed_version: 0,
                        violations,
                    });
                }
                for batch in schema::MIGRATIONS {
                    for ddl in *batch {
                        self.conn.execute(ddl).map_err(|error| LedgerError::Sql {
                            context: "initialize schema".to_string(),
                            detail: format!(
                                "{error} while executing: {}",
                                ddl.get(..60).unwrap_or(ddl)
                            ),
                        })?;
                    }
                }
                self.seed_instance_id_if_missing()?;
                let _ = self.read_current_instance_id()?;
                self.verify_evidence_semantic_bindings()?;
                self.backfill_and_verify_op_content_identities()?;
                self.backfill_and_verify_tune_content_identities()?;
                self.conn
                    .execute(&format!("PRAGMA user_version = {SCHEMA_VERSION}"))
                    .map_err(|error| sql_err("init: set user_version", &error))?;
                self.conn
                    .commit_transaction()
                    .map_err(|error| sql_err("init: commit", &error))
            })();
            if let Err(error) = init {
                let _ = self.conn.rollback_transaction();
                return Err(error);
            }
            return Ok(());
        }
        self.attest_schema(found)?;
        if found >= 4 {
            // Identity became mandatory in v4. Validate it before any later
            // schema marker advances, so a corrupt v4 file remains v4.
            let _ = self.read_current_instance_id()?;
        }
        let start = usize::try_from(found).unwrap_or(usize::MAX);
        for (step, batch) in schema::MIGRATIONS.iter().enumerate().skip(start) {
            let target = step + 1;
            self.conn
                .begin_transaction()
                .map_err(|e| sql_err("migrate: begin", &e))?;
            let migration = (|| {
                if target == 5 {
                    // Repeat the v4 identity check inside the transaction that
                    // installs the guards, closing the check-to-trigger window.
                    let _ = self.read_current_instance_id()?;
                }
                for ddl in *batch {
                    let already_applied = schema::RECOVERABLE_ADDED_COLUMNS
                        .iter()
                        .find(|column| column.ddl == *ddl)
                        .map_or(Ok(false), |column| {
                            self.recoverable_column_is_present(target, column)
                        })?;
                    if already_applied {
                        continue;
                    }
                    self.conn.execute(ddl).map_err(|error| LedgerError::Sql {
                        context: format!("migrate to v{target}"),
                        detail: format!(
                            "{error} while executing: {}",
                            ddl.get(..60).unwrap_or(ddl)
                        ),
                    })?;
                }
                if target == 4 {
                    self.seed_instance_id_if_missing()?;
                    let _ = self.read_current_instance_id()?;
                }
                if target == 8 {
                    // V8 backfills a redundant immutable discovery witness from
                    // v6/v7 claims. Authenticate every source claim before the
                    // version marker commits so valid-looking pre-migration
                    // semantic corruption cannot be copied into both indexes.
                    self.verify_session_claim_discovery_backfill()?;
                }
                if target == 14 {
                    // V14 copies only the exact raw artifact hash into the
                    // typed ContentId sidecar. Re-hash every retained source
                    // and verify the complete one-to-one mapping before the
                    // marker commits; failure rolls back DDL and data alike.
                    self.verify_artifact_content_identity_backfill()?;
                }
                if target == 15 {
                    // V15 carries the already-authenticated artifact ContentId
                    // through every role-qualified lineage edge. Audit both
                    // layers before the marker commits.
                    self.verify_edge_content_identity_backfill()?;
                }
                if target == 16 {
                    // V16 infers no bindings. Authenticate any exact rows
                    // present ahead of a stale marker and every prerequisite
                    // identity layer before committing the marker.
                    self.verify_artifact_semantic_bindings()?;
                }
                if target == 17 {
                    // V17 infers no evidence semantics. Authenticate exact
                    // preapplied bindings and every prerequisite layer before
                    // committing a stale marker.
                    self.verify_evidence_semantic_bindings()?;
                }
                if target == 18 {
                    // V18 computes no semantic meaning. Hash each bounded
                    // frozen operation field in Rust, retain the raw-content
                    // types separately, and verify every source/sidecar pair
                    // before the marker commits.
                    self.backfill_and_verify_op_content_identities()?;
                }
                if target == 19 {
                    // V19 hashes only the exact bounded cache key/value bytes.
                    // Re-hash every field independently and refuse a partial
                    // or semantically guessed cache migration.
                    self.backfill_and_verify_tune_content_identities()?;
                }
                self.conn
                    .execute(&format!("PRAGMA user_version = {target}"))
                    .map_err(|error| sql_err("migrate: set user_version", &error))?;
                self.conn
                    .commit_transaction()
                    .map_err(|error| sql_err("migrate: commit", &error))
            })();
            if let Err(error) = migration {
                let _ = self.conn.rollback_transaction();
                return Err(error);
            }
        }
        Ok(())
    }

    /// Re-attest the exact identity mutation guards without rebuilding and
    /// diffing the whole schema. The lookup is bounded to the three reserved
    /// names plus one duplicate sentinel row; healthy calls issue one small
    /// sqlite_master query and allocate no unbounded database-derived state.
    fn attest_identity_guards(&self) -> Result<(), LedgerError> {
        let rows = self
            .conn
            .query(
                "SELECT type, name, COALESCE(sql, '') FROM sqlite_master \
                 WHERE name IN ( \
                   'trg_ledger_identity_immutable_update', \
                   'trg_ledger_identity_immutable_delete', \
                   'trg_ledger_identity_immutable_reinsert' \
                 ) ORDER BY name LIMIT 4",
            )
            .map_err(|error| sql_err("attest ledger identity guards", &error))?;

        let mut actual = std::collections::BTreeMap::new();
        let mut violations = Vec::new();
        for row in &rows {
            let kind = row_text(row, 0, "attest identity guard type")?;
            let name = row_text(row, 1, "attest identity guard name")?;
            let sql = row_text(row, 2, "attest identity guard SQL")?;
            if actual
                .insert(name.clone(), (kind, normalize_schema_sql(&sql)))
                .is_some()
            {
                violations.push(format!("duplicate schema object `{name}`"));
            }
        }

        for (name, ddl) in LEDGER_IDENTITY_GUARDS {
            let expected_sql = normalize_schema_sql(ddl);
            match actual.get(name) {
                None => violations.push(format!("missing trigger `{name}`")),
                Some((kind, _)) if kind.as_str() != "trigger" => {
                    violations.push(format!(
                        "identity guard `{name}` is a {kind}, not a trigger"
                    ));
                }
                Some((_, sql)) if sql != &expected_sql => violations.push(format!(
                    "trigger `{name}` differs from the shipped identity-guard definition"
                )),
                Some(_) => {}
            }
        }

        if violations.is_empty() {
            return Ok(());
        }
        Err(LedgerError::SchemaMismatch {
            claimed_version: self.schema_version()?,
            violations,
        })
    }

    fn seed_instance_id_if_missing(&self) -> Result<(), LedgerError> {
        let existing = self
            .conn
            .query("SELECT singleton, instance_id FROM ledger_identity LIMIT 2")
            .map_err(|error| sql_err("inspect ledger instance identity", &error))?;
        if !existing.is_empty() {
            let _ = self.read_current_instance_id()?;
            return Ok(());
        }
        let instance_id = fresh_ledger_instance_id()?;
        self.conn
            .prepare(
                "INSERT INTO ledger_identity(singleton, instance_id) \
                 SELECT 1, ?1 WHERE NOT EXISTS \
                 (SELECT 1 FROM ledger_identity WHERE singleton = 1)",
            )
            .map_err(|error| sql_err("seed ledger instance identity: prepare", &error))?
            .execute_with_params(&[blob_param(&instance_id.0)])
            .map_err(|error| sql_err("seed ledger instance identity", &error))?;
        Ok(())
    }

    fn read_current_instance_id(&self) -> Result<LedgerInstanceId, LedgerError> {
        let rows = self
            .conn
            .query("SELECT singleton, instance_id FROM ledger_identity LIMIT 2")
            .map_err(|error| sql_err("read ledger instance identity", &error))?;
        decode_ledger_instance_id(
            rows.len(),
            rows.first().and_then(|row| row.get(0)),
            rows.first().and_then(|row| row.get(1)),
        )
    }

    fn recoverable_column_is_present(
        &self,
        target: usize,
        expected: &schema::RecoverableAddedColumn,
    ) -> Result<bool, LedgerError> {
        let rows = self
            .conn
            .query(&format!("PRAGMA table_info({})", expected.table))
            .map_err(|error| sql_err("migrate: inspect recoverable column", &error))?;
        for row in &rows {
            let name = match row.get(1) {
                Some(SqliteValue::Text(value)) => value.as_str(),
                other => {
                    return Err(LedgerError::Sql {
                        context: format!("migrate to v{target}"),
                        detail: format!(
                            "PRAGMA table_info({}) returned a non-TEXT column name: {other:?}",
                            expected.table
                        ),
                    });
                }
            };
            if name != expected.name {
                continue;
            }
            let declared_type = match row.get(2) {
                Some(SqliteValue::Text(value)) => value.as_str(),
                other => {
                    return Err(LedgerError::Sql {
                        context: format!("migrate to v{target}"),
                        detail: format!(
                            "recoverable column {}.{} has non-TEXT declared type {other:?}",
                            expected.table, expected.name
                        ),
                    });
                }
            };
            let not_null = row_i64(row, 3, "migrate: inspect column not-null")? != 0;
            let default_sql = match row.get(4) {
                Some(SqliteValue::Null) => None,
                Some(SqliteValue::Text(value)) => Some(value.as_str()),
                other => {
                    return Err(LedgerError::Sql {
                        context: format!("migrate to v{target}"),
                        detail: format!(
                            "recoverable column {}.{} has invalid default metadata {other:?}",
                            expected.table, expected.name
                        ),
                    });
                }
            };
            let primary_key = row_i64(row, 5, "migrate: inspect column primary-key")? != 0;
            if declared_type.eq_ignore_ascii_case(expected.declared_type)
                && not_null == expected.not_null
                && default_sql == expected.default_sql
                && primary_key == expected.primary_key
            {
                return Ok(true);
            }
            return Err(LedgerError::Sql {
                context: format!("migrate to v{target}"),
                detail: format!(
                    "existing column {}.{} does not match the recoverable migration definition: \
                     type={declared_type:?}, not_null={not_null}, default={default_sql:?}, \
                     primary_key={primary_key}",
                    expected.table, expected.name
                ),
            });
        }
        Ok(false)
    }

    /// ATTEST the on-disk objects against a REFERENCE database built by
    /// replaying the shipped DDL up to `claimed` (bead gp3.18). Compares
    /// sqlite_master object-for-object (tables, indexes, triggers, views —
    /// the stored SQL text covers foreign keys, CHECKs, and STRICT) and
    /// every table column-for-column via PRAGMA table_info (name, declared
    /// type/affinity, not-null, default, primary key). RECOVERY TOLERANCE
    /// (the tt_001b crash-window contract, generalizing
    /// RECOVERABLE_ADDED_COLUMNS): objects or columns beyond `claimed` are
    /// tolerated IFF they match the CURRENT shipped schema exactly — a
    /// committed-DDL-but-stale-marker file heals, while any divergent
    /// early object fails closed.
    fn attest_schema(&self, claimed: i64) -> Result<(), LedgerError> {
        let steps = usize::try_from(claimed).unwrap_or(usize::MAX);
        let at_claimed = reference_connection(steps)?;
        let at_current = reference_connection(schema::MIGRATIONS.len())?;
        let mut violations = Vec::new();
        let expected = schema_objects(&at_claimed)?;
        let current = schema_objects(&at_current)?;
        let actual = schema_objects(&self.conn)?;
        for (key, sql) in &expected {
            match actual.get(key) {
                None => violations.push(format!("missing {} `{}`", key.0, key.1)),
                Some(found) if found != sql && Some(found) != current.get(key) => violations.push(
                    format!("{} `{}` differs from the shipped definition", key.0, key.1),
                ),
                Some(_) => {}
            }
        }
        for (key, sql) in &actual {
            if !expected.contains_key(key) && current.get(key) != Some(sql) {
                violations.push(format!(
                    "unexpected {} `{}` (conflicting or foreign object)",
                    key.0, key.1
                ));
            }
        }
        // Column-level attestation for every expected table: COLUMN-level
        // diagnostics, with the same future-form tolerance (a column from a
        // committed-but-unmarked later batch must match its shipped
        // definition exactly).
        for (kind, table) in expected.keys() {
            if kind != "table" || !actual.contains_key(&(kind.clone(), table.clone())) {
                continue;
            }
            let want = table_columns(&at_claimed, table)?;
            let full = table_columns(&at_current, table)?;
            let have = table_columns(&self.conn, table)?;
            for col in &want {
                if !have.contains(col) {
                    violations.push(format!("table `{table}`: missing or altered column {col}"));
                }
            }
            for col in &have {
                if !want.contains(col) && !full.contains(col) {
                    violations.push(format!("table `{table}`: unexpected column {col}"));
                }
            }
        }
        if violations.is_empty() {
            Ok(())
        } else {
            violations.sort();
            Err(LedgerError::SchemaMismatch {
                claimed_version: claimed,
                violations,
            })
        }
    }

    /// `json_valid` check through the same engine that enforces the schema
    /// CHECKs, so pre-validation and enforcement can never disagree.
    fn json_valid(&self, s: &str) -> Result<bool, LedgerError> {
        let rows = self
            .conn
            .query_with_params("SELECT json_valid(?1)", &[text_param(s)])
            .map_err(|e| sql_err("json_valid", &e))?;
        let row = rows.first().ok_or_else(|| LedgerError::Sql {
            context: "json_valid".to_string(),
            detail: "no row returned".to_string(),
        })?;
        Ok(row_i64(row, 0, "json_valid")? == 1)
    }

    fn require_json(&self, field: &str, value: &str, explicit: bool) -> Result<(), LedgerError> {
        let problem = if value.trim().is_empty() {
            Some("empty string; supply a JSON value".to_string())
        } else if !self.json_valid(value)? {
            Some(format!(
                "not valid JSON: {:?}",
                value.get(..40).unwrap_or(value)
            ))
        } else {
            None
        };
        match problem {
            None => Ok(()),
            Some(problem) if explicit => Err(LedgerError::MissingExplicit {
                field: field.to_string(),
                problem,
            }),
            Some(problem) => Err(LedgerError::Invalid {
                field: field.to_string(),
                problem,
            }),
        }
    }

    #[allow(clippy::unused_self)] // Shared impl modules call this stateless validator as a Ledger method.
    fn require_op_field_bound(
        &self,
        field: &str,
        len: usize,
        max_bytes: usize,
        explicit: bool,
    ) -> Result<(), LedgerError> {
        if len <= max_bytes {
            return Ok(());
        }
        let problem = format!(
            "{len} bytes exceeds the {max_bytes}-byte operation-field limit; shorten the frozen provenance value"
        );
        if explicit {
            Err(LedgerError::MissingExplicit {
                field: field.to_string(),
                problem,
            })
        } else {
            Err(LedgerError::Invalid {
                field: field.to_string(),
                problem,
            })
        }
    }

    fn require_bounded_op_json(
        &self,
        field: &str,
        value: &str,
        max_bytes: usize,
        explicit: bool,
    ) -> Result<(), LedgerError> {
        self.require_op_field_bound(field, value.len(), max_bytes, explicit)?;
        self.require_json(field, value, explicit)
    }

    fn require_tune_json(
        &self,
        field: &str,
        value: &str,
        max_bytes: usize,
    ) -> Result<(), LedgerError> {
        if value.len() > max_bytes {
            return Err(LedgerError::Invalid {
                field: field.to_string(),
                problem: format!(
                    "{} bytes exceeds the {max_bytes}-byte tune JSON limit",
                    value.len()
                ),
            });
        }
        self.require_json(field, value, false)
    }

    fn require_tune_row(
        &self,
        kernel: &str,
        shape_class: &str,
        machine: &[u8],
        params: &str,
        measured: &str,
    ) -> Result<(), LedgerError> {
        validate_tune_identity("kernel", kernel, MAX_TUNE_KERNEL_BYTES)?;
        validate_tune_identity("shape_class", shape_class, MAX_TUNE_SHAPE_CLASS_BYTES)?;
        validate_tune_machine(machine)?;
        self.require_tune_json("params", params, MAX_TUNE_PARAMS_BYTES)?;
        self.require_tune_json("measured", measured, MAX_TUNE_MEASURED_BYTES)
    }

    // -- transactions -------------------------------------------------------

    /// Begin an explicit transaction (for atomic op+edges+metrics groups).
    ///
    /// # Errors
    /// [`LedgerError::Busy`] on contention; [`LedgerError::Sql`] otherwise.
    pub fn begin(&self) -> Result<(), LedgerError> {
        self.conn
            .begin_transaction()
            .map_err(|e| sql_err("begin", &e))
    }

    /// Commit the explicit transaction.
    ///
    /// # Errors
    /// [`LedgerError::Busy`] on commit-time conflict (retry the whole group).
    pub fn commit(&self) -> Result<(), LedgerError> {
        self.conn
            .commit_transaction()
            .map_err(|e| sql_err("commit", &e))
    }

    /// Roll back the explicit transaction.
    ///
    /// # Errors
    /// [`LedgerError::Sql`] on engine failure.
    pub fn rollback(&self) -> Result<(), LedgerError> {
        self.conn
            .rollback_transaction()
            .map_err(|e| sql_err("rollback", &e))
    }

    /// True while an explicit transaction is open.
    #[must_use]
    pub fn in_transaction(&self) -> bool {
        self.conn.in_transaction()
    }

    // -- artifacts ----------------------------------------------------------

    /// Store `bytes` content-addressed. Identical bytes dedupe to one row,
    /// across forks, forever; the receipt says whether this call created a
    /// new row. Artifacts larger than [`STORAGE_CHUNK_LEN`] are stored as
    /// chunk rows (fsqlite has no incremental-blob API; RAM per row stays
    /// bounded).
    ///
    /// # Errors
    /// [`LedgerError::Invalid`] for empty `kind`/malformed `meta`;
    /// [`LedgerError::Busy`]/[`LedgerError::Sql`] from the engine.
    pub fn put_artifact(
        &self,
        kind: &str,
        bytes: &[u8],
        meta: Option<&str>,
    ) -> Result<PutReceipt, LedgerError> {
        self.validate_artifact_inputs(kind, meta)?;
        let h = hash_bytes(bytes);
        if let Some(info) = self.artifact_info(&h)? {
            self.verify_dedupe_candidate(&h, &info, kind, meta)?;
            return Ok(PutReceipt {
                hash: h,
                len: info.len,
                deduped: true,
                chunked: info.chunk_count > 0,
            });
        }
        let len = bytes.len() as u64;
        if bytes.len() <= STORAGE_CHUNK_LEN {
            return self.insert_inline_artifact(&h, kind, bytes, meta);
        }
        // Chunked path, atomically.
        let owns_txn = !self.conn.in_transaction();
        if owns_txn {
            self.begin()?;
        }
        let result = self.insert_chunked_artifact(&h, kind, bytes, meta);
        match (&result, owns_txn) {
            (Ok(_), true) => {
                if let Err(e) = self.commit() {
                    let _ = self.rollback();
                    return Err(e);
                }
            }
            (Err(_), true) => {
                let _ = self.rollback();
            }
            _ => {}
        }
        match result? {
            ChunkedArtifactInsert::Inserted => Ok(PutReceipt {
                hash: h,
                len,
                deduped: false,
                chunked: true,
            }),
            ChunkedArtifactInsert::Deduped(info) => Ok(PutReceipt {
                hash: h,
                len: info.len,
                deduped: true,
                chunked: info.chunk_count > 0,
            }),
        }
    }

    fn validate_artifact_inputs(&self, kind: &str, meta: Option<&str>) -> Result<(), LedgerError> {
        if kind.is_empty() {
            return Err(LedgerError::Invalid {
                field: "kind".to_string(),
                problem: "empty; name the artifact kind (e.g. \"field\", \"mesh\")".to_string(),
            });
        }
        if kind.len() > MAX_ARTIFACT_KIND_BYTES {
            return Err(LedgerError::Invalid {
                field: "kind".to_string(),
                problem: format!(
                    "UTF-8 encoding exceeds the {MAX_ARTIFACT_KIND_BYTES}-byte artifact-kind limit"
                ),
            });
        }
        if let Some(m) = meta {
            if m.len() > MAX_ARTIFACT_META_BYTES {
                return Err(LedgerError::Invalid {
                    field: "meta".to_string(),
                    problem: format!(
                        "UTF-8 encoding exceeds the {MAX_ARTIFACT_META_BYTES}-byte artifact-metadata limit"
                    ),
                });
            }
            self.require_json("meta", m, false)?;
        }
        Ok(())
    }

    /// Envelope agreement gate at every dedupe site (bead gp3.19): a
    /// byte-identical artifact dedupes only when the offered `kind`
    /// matches exactly and the offered metadata (when a claim is made)
    /// canonically equals the stored metadata. Offering `meta: None`
    /// makes NO claim and accepts the stored envelope (the streaming
    /// dedupe contract); offering metadata against a row stored without
    /// any is a conflict — silent claim-dropping is the bug this gate
    /// removes.
    fn attest_artifact_envelope(
        &self,
        h: &ContentHash,
        info: &ArtifactInfo,
        kind: &str,
        meta: Option<&str>,
    ) -> Result<(), LedgerError> {
        if kind != info.kind {
            return Err(LedgerError::ArtifactEnvelopeConflict {
                hash_hex: h.to_hex(),
                field: "kind",
                stored: info.kind.clone(),
                offered: kind.to_string(),
            });
        }
        let Some(offered) = meta else {
            return Ok(());
        };
        let conflict = |stored: String| LedgerError::ArtifactEnvelopeConflict {
            hash_hex: h.to_hex(),
            field: "meta",
            stored,
            offered: offered.to_string(),
        };
        let Some(stored) = info.meta.as_deref() else {
            return Err(conflict("<none>".to_string()));
        };
        // Canonical comparison through the engine (whitespace-insensitive;
        // key order remains significant, as documented in CONTRACT.md).
        let rows = self
            .conn
            .query_with_params(
                "SELECT json(?1) = json(?2)",
                &[text_param(stored), text_param(offered)],
            )
            .map_err(|e| sql_err("artifact meta compare", &e))?;
        let equal = rows
            .first()
            .map_or(Ok(0), |row| row_i64(row, 0, "artifact meta compare"))?;
        if equal == 1 {
            Ok(())
        } else {
            Err(conflict(stored.to_string()))
        }
    }

    /// Refuse dedupe unless the existing storage still hashes to its key.
    /// Envelope agreement alone cannot make corrupted bytes trustworthy.
    fn verify_dedupe_candidate(
        &self,
        h: &ContentHash,
        info: &ArtifactInfo,
        kind: &str,
        meta: Option<&str>,
    ) -> Result<(), LedgerError> {
        match self.read_artifact_chunks(h, &mut |_| {})? {
            Some(streamed) if streamed == info.len => {}
            Some(streamed) => {
                return Err(LedgerError::Corrupt {
                    hash_hex: h.to_hex(),
                    detail: format!(
                        "dedupe verification streamed {streamed} bytes but metadata records {}",
                        info.len
                    ),
                });
            }
            None => {
                return Err(LedgerError::Corrupt {
                    hash_hex: h.to_hex(),
                    detail: "dedupe candidate disappeared before its bytes could be verified"
                        .to_string(),
                });
            }
        }
        self.attest_artifact_envelope(h, info, kind, meta)
    }

    fn insert_inline_artifact(
        &self,
        h: &ContentHash,
        kind: &str,
        bytes: &[u8],
        meta: Option<&str>,
    ) -> Result<PutReceipt, LedgerError> {
        let insert = self
            .conn
            .prepare(
                "INSERT INTO artifacts(hash, kind, bytes, len, chunk_count, meta, created_at) \
                 VALUES (?1, ?2, ?3, ?4, 0, ?5, ?6)",
            )
            .map_err(|e| sql_err("artifact insert prepare", &e))?
            .execute_with_params(&[
                blob_param(h.as_bytes()),
                text_param(kind),
                blob_param(bytes),
                SqliteValue::Integer(int_from_usize(bytes.len())),
                opt_text_param(meta),
                SqliteValue::Integer(now_wall_ns()),
            ]);
        match insert {
            Ok(_) => Ok(PutReceipt {
                hash: *h,
                len: bytes.len() as u64,
                deduped: false,
                chunked: false,
            }),
            // A concurrent writer stored the same content first: that IS the
            // dedupe contract, not an error — but only under an AGREEING
            // envelope (gp3.19), so re-read theirs and attest.
            Err(e) if is_duplicate_key(&e) => {
                let info = self.artifact_info(h)?.ok_or_else(|| LedgerError::Corrupt {
                    hash_hex: h.to_hex(),
                    detail: "duplicate artifact key exists without readable artifact metadata"
                        .to_string(),
                })?;
                self.verify_dedupe_candidate(h, &info, kind, meta)?;
                Ok(PutReceipt {
                    hash: *h,
                    len: info.len,
                    deduped: true,
                    chunked: info.chunk_count > 0,
                })
            }
            Err(e) => Err(sql_err("artifact insert", &e)),
        }
    }

    fn insert_chunked_artifact(
        &self,
        h: &ContentHash,
        kind: &str,
        bytes: &[u8],
        meta: Option<&str>,
    ) -> Result<ChunkedArtifactInsert, LedgerError> {
        let chunks = bytes.chunks(STORAGE_CHUNK_LEN);
        let chunk_count = int_from_usize(chunks.len());
        for (seq, chunk) in chunks.enumerate() {
            let insert = self
                .conn
                .prepare("INSERT INTO artifact_chunks(hash, seq, bytes) VALUES (?1, ?2, ?3)")
                .map_err(|e| sql_err("chunk insert prepare", &e))?
                .execute_with_params(&[
                    blob_param(h.as_bytes()),
                    SqliteValue::Integer(int_from_usize(seq)),
                    blob_param(chunk),
                ]);
            match insert {
                Ok(_) => {}
                Err(e) if is_duplicate_key(&e) => {
                    // Concurrent identical store; the other writer wins —
                    // if their committed bytes and envelope agree. Restore
                    // any holes this attempt filled before evaluating the
                    // pre-existing candidate.
                    if seq > 0 {
                        self.conn
                            .prepare("DELETE FROM artifact_chunks WHERE hash = ?1 AND seq < ?2")
                            .map_err(|error| sql_err("chunk prefix cleanup prepare", &error))?
                            .execute_with_params(&[
                                blob_param(h.as_bytes()),
                                SqliteValue::Integer(int_from_usize(seq)),
                            ])
                            .map_err(|error| sql_err("chunk prefix cleanup", &error))?;
                    }
                    let info = self.artifact_info(h)?.ok_or_else(|| LedgerError::Corrupt {
                        hash_hex: h.to_hex(),
                        detail: "duplicate chunk key exists without artifact metadata".to_string(),
                    })?;
                    self.verify_dedupe_candidate(h, &info, kind, meta)?;
                    return Ok(ChunkedArtifactInsert::Deduped(info));
                }
                Err(e) => return Err(sql_err("chunk insert", &e)),
            }
        }
        self.conn
            .prepare(
                "INSERT INTO artifacts(hash, kind, bytes, len, chunk_count, meta, created_at) \
                 VALUES (?1, ?2, NULL, ?3, ?4, ?5, ?6)",
            )
            .map_err(|e| sql_err("artifact insert prepare", &e))?
            .execute_with_params(&[
                blob_param(h.as_bytes()),
                text_param(kind),
                SqliteValue::Integer(int_from_usize(bytes.len())),
                SqliteValue::Integer(chunk_count),
                opt_text_param(meta),
                SqliteValue::Integer(now_wall_ns()),
            ])
            .map_err(|e| sql_err("artifact insert", &e))?;
        Ok(ChunkedArtifactInsert::Inserted)
    }

    /// Begin a streaming artifact write (for fields too large to hold in
    /// memory). The writer owns a transaction: a crash or drop before
    /// `finish` leaves zero residue.
    ///
    /// # Errors
    /// [`LedgerError::WriterInTransaction`] if an explicit transaction is
    /// open; engine errors otherwise.
    pub fn artifact_writer(&self, kind: &str) -> Result<ArtifactWriter<'_>, LedgerError> {
        self.validate_artifact_inputs(kind, None)?;
        if self.conn.in_transaction() {
            return Err(LedgerError::WriterInTransaction);
        }
        self.begin()?;
        Ok(ArtifactWriter {
            ledger: self,
            kind: kind.to_string(),
            hasher: Blake3::new(),
            provisional: provisional_key(),
            next_seq: 0,
            buf: Vec::new(),
            len: 0,
            finished: false,
        })
    }

    fn artifact_envelope_is_bounded(&self, h: &ContentHash) -> Result<bool, LedgerError> {
        let preflight = self
            .conn
            .query_with_params(
                "SELECT typeof(kind), length(CAST(kind AS BLOB)), typeof(meta), \
                 CASE WHEN meta IS NULL THEN 0 ELSE length(CAST(meta AS BLOB)) END \
                 FROM artifacts WHERE hash = ?1",
                &[blob_param(h.as_bytes())],
            )
            .map_err(|e| sql_err("artifact_info envelope preflight", &e))?;
        let Some(preflight) = preflight.first() else {
            return Ok(false);
        };
        let envelope_detail = (|| {
            let kind_type = match preflight.get(0) {
                Some(SqliteValue::Text(value)) => value.as_str(),
                other => {
                    return Err(format!(
                        "artifact kind type preflight expected TEXT, got {other:?}"
                    ));
                }
            };
            let meta_type = match preflight.get(2) {
                Some(SqliteValue::Text(value)) => value.as_str(),
                other => {
                    return Err(format!(
                        "artifact metadata type preflight expected TEXT, got {other:?}"
                    ));
                }
            };
            let kind_len = storage_u64(preflight, 1, "artifact kind byte length")?;
            let meta_len = storage_u64(preflight, 3, "artifact metadata byte length")?;
            if kind_type != "text" {
                return Err(format!("artifact kind must be TEXT, found {kind_type}"));
            }
            if kind_len == 0 || kind_len > MAX_ARTIFACT_KIND_BYTES as u64 {
                return Err(format!(
                    "artifact kind byte length {kind_len} is outside 1..={MAX_ARTIFACT_KIND_BYTES}"
                ));
            }
            if meta_type != "null" && meta_type != "text" {
                return Err(format!(
                    "artifact metadata must be TEXT or NULL, found {meta_type}"
                ));
            }
            if meta_len > MAX_ARTIFACT_META_BYTES as u64 {
                return Err(format!(
                    "artifact metadata byte length {meta_len} exceeds {MAX_ARTIFACT_META_BYTES}"
                ));
            }
            Ok(())
        })();
        envelope_detail.map_err(|detail| LedgerError::Corrupt {
            hash_hex: h.to_hex(),
            detail,
        })?;
        Ok(true)
    }

    /// Metadata for one artifact, if present.
    ///
    /// # Errors
    /// Engine errors or [`LedgerError::Corrupt`] when the stored envelope is
    /// malformed or exceeds its materialization bounds; absence is `Ok(None)`.
    pub fn artifact_info(&self, h: &ContentHash) -> Result<Option<ArtifactInfo>, LedgerError> {
        if !self.artifact_envelope_is_bounded(h)? {
            return Ok(None);
        }
        let rows = self
            .conn
            .query_with_params(
                "SELECT kind, len, chunk_count, meta, created_at FROM artifacts \
                 WHERE hash = ?1 AND typeof(kind) = 'text' \
                 AND length(CAST(kind AS BLOB)) BETWEEN 1 AND ?2 \
                 AND (meta IS NULL OR (typeof(meta) = 'text' \
                      AND length(CAST(meta AS BLOB)) <= ?3))",
                &[
                    blob_param(h.as_bytes()),
                    SqliteValue::Integer(int_from_usize(MAX_ARTIFACT_KIND_BYTES)),
                    SqliteValue::Integer(int_from_usize(MAX_ARTIFACT_META_BYTES)),
                ],
            )
            .map_err(|e| sql_err("artifact_info", &e))?;
        let row = rows.first().ok_or_else(|| LedgerError::Corrupt {
            hash_hex: h.to_hex(),
            detail: "artifact envelope disappeared or exceeded its bound after preflight"
                .to_string(),
        })?;
        let kind = match row.get(0) {
            Some(SqliteValue::Text(t)) => t.as_str().to_string(),
            other => {
                return Err(LedgerError::Sql {
                    context: "artifact_info".to_string(),
                    detail: format!("kind: expected TEXT, got {other:?}"),
                });
            }
        };
        let meta = match row.get(3) {
            Some(SqliteValue::Text(t)) => Some(t.as_str().to_string()),
            Some(SqliteValue::Null) | None => None,
            other => {
                return Err(LedgerError::Sql {
                    context: "artifact_info".to_string(),
                    detail: format!("meta: expected TEXT or NULL, got {other:?}"),
                });
            }
        };
        if let Some(meta) = meta.as_deref() {
            match self.require_json("artifact.meta", meta, false) {
                Ok(()) => {}
                Err(LedgerError::Invalid { .. }) => {
                    return Err(LedgerError::Corrupt {
                        hash_hex: h.to_hex(),
                        detail: "artifact metadata is not valid JSON".to_string(),
                    });
                }
                Err(error) => return Err(error),
            }
        }
        Ok(Some(ArtifactInfo {
            hash: *h,
            kind,
            len: row_u64(row, 1, "artifact_info.len")?,
            chunk_count: row_u64(row, 2, "artifact_info.chunk_count")?,
            meta,
            created_at: row_i64(row, 4, "artifact_info.created_at")?,
        }))
    }

    fn preflight_inline_artifact(
        &self,
        h: &ContentHash,
        info: &ArtifactInfo,
    ) -> Result<(), LedgerError> {
        let rows = self
            .conn
            .query_with_params(
                "SELECT CASE WHEN bytes IS NULL THEN 0 ELSE 1 END, \
                 COALESCE(length(bytes), -1) FROM artifacts WHERE hash = ?1",
                &[blob_param(h.as_bytes())],
            )
            .map_err(|e| sql_err("inline artifact preflight", &e))?;
        let row = rows.first().ok_or_else(|| LedgerError::Corrupt {
            hash_hex: h.to_hex(),
            detail: "artifact metadata disappeared before storage preflight".to_string(),
        })?;
        let detail = (|| {
            let present = storage_i64(row, 0, "inline BLOB presence")?;
            let actual_len = storage_i64(row, 1, "inline BLOB length")?;
            if present != 1 || actual_len < 0 {
                return Err("inline artifact has no BLOB value".to_string());
            }
            let actual_len = u64::try_from(actual_len)
                .map_err(|_| "inline BLOB length is outside the ledger domain".to_string())?;
            let bound = STORAGE_CHUNK_LEN as u64;
            if info.len > bound || actual_len > bound {
                return Err(format!(
                    "inline value exceeds the {STORAGE_CHUNK_LEN}-byte storage bound: recorded {}, found {actual_len}",
                    info.len
                ));
            }
            if actual_len != info.len {
                return Err(format!(
                    "inline length mismatch: recorded {}, found {actual_len}",
                    info.len
                ));
            }
            Ok(())
        })();
        detail.map_err(|detail| LedgerError::Corrupt {
            hash_hex: h.to_hex(),
            detail,
        })
    }

    fn preflight_chunked_artifact(
        &self,
        h: &ContentHash,
        info: &ArtifactInfo,
    ) -> Result<(), LedgerError> {
        let row = self
            .conn
            .query_row_with_params(
                "SELECT COUNT(*), COUNT(bytes), COALESCE(MIN(seq), -1), \
                 COALESCE(MAX(seq), -1), COALESCE(SUM(length(bytes)), 0), \
                 COALESCE(MAX(length(bytes)), 0) \
                 FROM artifact_chunks WHERE hash = ?1",
                &[blob_param(h.as_bytes())],
            )
            .map_err(|e| sql_err("chunked artifact preflight", &e))?;
        ArtifactChunkPreflight::from_row(&row)
            .and_then(|preflight| preflight.validate(info))
            .map_err(|detail| LedgerError::Corrupt {
                hash_hex: h.to_hex(),
                detail,
            })
    }

    /// Fetch an artifact's full bytes (assembles chunked storage in memory —
    /// prefer [`Ledger::read_artifact_chunks`] for very large fields). The
    /// materializer reserves fallibly from bytes actually read; recorded
    /// length metadata is never used as an allocation request.
    ///
    /// # Errors
    /// Engine errors, [`LedgerError::Corrupt`] when storage does not match its
    /// recorded shape, or [`LedgerError::Invalid`] when memory cannot be
    /// reserved for the result; absence is `Ok(None)`.
    pub fn get_artifact(&self, h: &ContentHash) -> Result<Option<Vec<u8>>, LedgerError> {
        self.note_read_query();
        let Some(info) = self.artifact_info(h)? else {
            return Ok(None);
        };
        self.materialize_artifact_with_info(h, &info).map(Some)
    }

    /// Fetch an artifact only when its metadata-declared length is within the
    /// caller's explicit materialization budget. The length comparison occurs
    /// before any payload callback or result-buffer allocation.
    ///
    /// # Errors
    /// The same errors as [`Ledger::get_artifact`], plus
    /// [`LedgerError::ArtifactReadLimit`] when stored metadata declares an
    /// artifact larger than `max_bytes`.
    pub fn get_artifact_bounded(
        &self,
        h: &ContentHash,
        max_bytes: u64,
    ) -> Result<Option<Vec<u8>>, LedgerError> {
        self.note_read_query();
        let Some(info) = self.artifact_info(h)? else {
            return Ok(None);
        };
        self.require_artifact_read_limit(h, &info, max_bytes)?;
        self.materialize_artifact_with_info(h, &info).map(Some)
    }

    #[allow(clippy::unused_self)] // Mirrors the method-shaped streaming/materialization validators.
    fn require_artifact_read_limit(
        &self,
        h: &ContentHash,
        info: &ArtifactInfo,
        max_bytes: u64,
    ) -> Result<(), LedgerError> {
        if info.len <= max_bytes {
            Ok(())
        } else {
            Err(LedgerError::ArtifactReadLimit {
                hash_hex: h.to_hex(),
                limit: max_bytes,
                observed: info.len,
            })
        }
    }

    fn materialize_artifact_with_info(
        &self,
        h: &ContentHash,
        info: &ArtifactInfo,
    ) -> Result<Vec<u8>, LedgerError> {
        let mut out = Vec::new();
        let mut allocation_error = None;
        self.read_artifact_chunks_with_info(h, info, &mut |chunk| {
            if allocation_error.is_some() {
                return;
            }
            if let Err(error) = out.try_reserve(chunk.len()) {
                allocation_error = Some(LedgerError::Invalid {
                    field: "artifact_bytes".to_string(),
                    problem: format!(
                        "could not reserve {} additional bytes ({error}); use \
                         Ledger::read_artifact_chunks to process the artifact incrementally",
                        chunk.len()
                    ),
                });
                return;
            }
            // `try_reserve` above established sufficient capacity, so this
            // copy cannot trigger another allocation.
            out.extend_from_slice(chunk);
        })?;
        if let Some(error) = allocation_error {
            return Err(error);
        }
        Ok(out)
    }

    /// Stream an artifact's bytes chunk-by-chunk without materializing the
    /// whole value. Returns the total length streamed, or `Ok(None)` if the
    /// artifact does not exist.
    ///
    /// Storage shape and row sizes are checked with metadata-only SQL before
    /// any BLOB is materialized. Each callback then receives only a prefix
    /// whose row sequence, count, and cumulative length are valid at that
    /// point. A concurrent mutation, a later row, or the final content-hash
    /// comparison can still return `Err`; callback side effects for an
    /// already-delivered prefix are not rolled back. In particular,
    /// same-length byte tampering is reported only after the full tampered
    /// prefix has been delivered and hashed. The first row that violates a
    /// recorded bound is not delivered. The callback has no error channel of
    /// its own, so callers that perform fallible work must record that failure
    /// and make subsequent invocations no-ops.
    ///
    /// # Errors
    /// Engine errors; [`LedgerError::Corrupt`] on malformed rows, non-dense
    /// sequences, oversized chunks, arithmetic overflow, disagreement with
    /// the recorded shape, or a content hash mismatch.
    pub fn read_artifact_chunks(
        &self,
        h: &ContentHash,
        f: &mut dyn FnMut(&[u8]),
    ) -> Result<Option<u64>, LedgerError> {
        let Some(info) = self.artifact_info(h)? else {
            return Ok(None);
        };
        self.read_artifact_chunks_with_info(h, &info, f).map(Some)
    }

    /// Stream an artifact only when its metadata-declared length is within the
    /// caller's explicit byte budget. A limit refusal occurs before storage
    /// preflight and before the callback can observe any payload byte.
    ///
    /// # Errors
    /// The same errors as [`Ledger::read_artifact_chunks`], plus
    /// [`LedgerError::ArtifactReadLimit`] when stored metadata declares an
    /// artifact larger than `max_bytes`.
    pub fn read_artifact_chunks_bounded(
        &self,
        h: &ContentHash,
        max_bytes: u64,
        f: &mut dyn FnMut(&[u8]),
    ) -> Result<Option<u64>, LedgerError> {
        let Some(info) = self.artifact_info(h)? else {
            return Ok(None);
        };
        self.require_artifact_read_limit(h, &info, max_bytes)?;
        self.read_artifact_chunks_with_info(h, &info, f).map(Some)
    }

    fn read_artifact_chunks_with_info(
        &self,
        h: &ContentHash,
        info: &ArtifactInfo,
        f: &mut dyn FnMut(&[u8]),
    ) -> Result<u64, LedgerError> {
        if info.chunk_count == 0 {
            self.preflight_inline_artifact(h, info)?;
            let rows = self
                .conn
                .query_with_params(
                    "SELECT bytes FROM artifacts WHERE hash = ?1 AND bytes IS NOT NULL \
                     AND length(bytes) <= ?2",
                    &[
                        blob_param(h.as_bytes()),
                        SqliteValue::Integer(int_from_usize(STORAGE_CHUNK_LEN)),
                    ],
                )
                .map_err(|e| sql_err("read_artifact_chunks", &e))?;
            let row = rows.first().ok_or_else(|| LedgerError::Corrupt {
                hash_hex: h.to_hex(),
                detail: "inline storage disappeared or exceeded its bound after preflight"
                    .to_string(),
            })?;
            let bytes =
                inline_artifact_bytes(row, info.len).map_err(|detail| LedgerError::Corrupt {
                    hash_hex: h.to_hex(),
                    detail,
                })?;
            let computed = hash_bytes(bytes);
            f(bytes);
            if computed != *h {
                return Err(LedgerError::Corrupt {
                    hash_hex: h.to_hex(),
                    detail: format!(
                        "content hash mismatch: computed {} from stored bytes",
                        computed.to_hex()
                    ),
                });
            }
            return Ok(info.len);
        }

        self.preflight_chunked_artifact(h, info)?;
        let mut validator = ArtifactChunkValidator::new(info);
        let mut hasher = Blake3::new();
        let mut corrupt_detail = None;
        self.conn
            .query_with_params_for_each(
                "SELECT seq, bytes FROM artifact_chunks WHERE hash = ?1 AND bytes IS NOT NULL \
                 AND length(bytes) <= ?2 ORDER BY seq",
                &[
                    blob_param(h.as_bytes()),
                    SqliteValue::Integer(int_from_usize(STORAGE_CHUNK_LEN)),
                ],
                |row| {
                    if corrupt_detail.is_some() {
                        return Ok(());
                    }
                    match validator.accept(row) {
                        Ok(bytes) => {
                            hasher.update(bytes);
                            f(bytes);
                        }
                        Err(detail) => corrupt_detail = Some(detail),
                    }
                    Ok(())
                },
            )
            .map_err(|e| sql_err("read_artifact_chunks", &e))?;
        if let Some(detail) = corrupt_detail {
            return Err(LedgerError::Corrupt {
                hash_hex: h.to_hex(),
                detail,
            });
        }
        let streamed = validator.finish().map_err(|detail| LedgerError::Corrupt {
            hash_hex: h.to_hex(),
            detail,
        })?;
        let computed = hasher.finalize();
        if computed != *h {
            return Err(LedgerError::Corrupt {
                hash_hex: h.to_hex(),
                detail: format!(
                    "content hash mismatch: computed {} from stored bytes",
                    computed.to_hex()
                ),
            });
        }
        Ok(streamed)
    }

    /// Re-hash every stored artifact against its recorded identity.
    /// Byte-level corruption fails LOUDLY here (Decalogue P9).
    ///
    /// # Errors
    /// Engine errors; corruption is reported in the result, not an `Err`.
    pub fn verify_artifact_integrity(&self) -> Result<IntegrityReport, LedgerError> {
        let rows = self
            .conn
            .query("SELECT hash FROM artifacts")
            .map_err(|e| sql_err("integrity scan", &e))?;
        let mut report = IntegrityReport::default();
        for row in &rows {
            let stored = match row.get(0) {
                Some(SqliteValue::Blob(b)) => ContentHash::from_slice(b),
                _ => None,
            };
            let Some(stored) = stored else {
                report.corrupted.push("<malformed hash column>".to_string());
                continue;
            };
            report.checked += 1;
            match self.read_artifact_chunks(&stored, &mut |_| {}) {
                Ok(Some(_)) => {}
                Ok(None) | Err(LedgerError::Corrupt { .. }) => {
                    report.corrupted.push(stored.to_hex());
                }
                Err(error) => return Err(error),
            }
        }
        Ok(report)
    }

    /// Test hook: overwrite one artifact's stored bytes so integrity checks
    /// can prove they fail loudly. Never call outside tests.
    ///
    /// # Errors
    /// Engine errors; [`LedgerError::NotFound`] if absent.
    pub fn corrupt_artifact_for_test(&self, h: &ContentHash) -> Result<(), LedgerError> {
        let Some(info) = self.artifact_info(h)? else {
            return Err(LedgerError::NotFound {
                what: format!("artifact {}", h.to_hex()),
            });
        };
        let (sql, params): (&str, Vec<SqliteValue>) = if info.chunk_count == 0 {
            (
                "UPDATE artifacts SET bytes = X'DEADBEEF', len = 4 WHERE hash = ?1",
                vec![blob_param(h.as_bytes())],
            )
        } else {
            (
                "UPDATE artifact_chunks SET bytes = X'DEADBEEF' WHERE hash = ?1 AND seq = 0",
                vec![blob_param(h.as_bytes())],
            )
        };
        self.conn
            .prepare(sql)
            .map_err(|e| sql_err("corrupt hook prepare", &e))?
            .execute_with_params(&params)
            .map_err(|e| sql_err("corrupt hook", &e))?;
        Ok(())
    }

    // -- ops and lineage ----------------------------------------------------

    /// Record the start of an op with its frozen Five Explicits (P4) on the
    /// main branch in deterministic mode (see [`Ledger::begin_op_on`] for
    /// forks and fast mode). The caller controls `t_start_ns` so
    /// deterministic replays can use logical time; [`now_wall_ns`] is the
    /// conventional wall-clock source.
    ///
    /// # Errors
    /// [`LedgerError::MissingExplicit`] naming the offending field;
    /// [`LedgerError::Invalid`] for an oversized session or malformed or
    /// oversized IR;
    /// engine errors otherwise.
    pub fn begin_op(
        &self,
        session: Option<&[u8]>,
        ir: &str,
        explicits: &FiveExplicits<'_>,
        t_start_ns: i64,
    ) -> Result<i64, LedgerError> {
        self.begin_op_on(
            MAIN_BRANCH,
            ExecMode::Deterministic,
            session,
            ir,
            explicits,
            t_start_ns,
        )
    }

    /// Record an op's outcome and seal its public-API lineage. Each op finishes
    /// exactly once; after this update commits, [`Ledger::link`] can no longer
    /// attach an input or output edge to the op.
    ///
    /// # Errors
    /// [`LedgerError::DoubleFinish`] on a second finish;
    /// [`LedgerError::NotFound`] for an unknown op;
    /// [`LedgerError::Invalid`] for an oversized or malformed diagnostic;
    /// [`LedgerError::OpCorrupt`] if a non-finishable stored row violates the
    /// bounded operation envelope.
    pub fn finish_op(
        &self,
        op: i64,
        outcome: OpOutcome,
        diag: Option<&str>,
        t_end_ns: i64,
    ) -> Result<(), LedgerError> {
        if let Some(d) = diag {
            self.require_bounded_op_json("diag", d, MAX_OP_DIAG_BYTES, false)?;
        }
        let affected = self
            .conn
            .prepare(
                "UPDATE ops SET t_end = ?1, outcome = ?2, diag = ?3 \
                 WHERE id = ?4 AND outcome IS NULL",
            )
            .map_err(|e| sql_err("op finish prepare", &e))?
            .execute_with_params(&[
                SqliteValue::Integer(t_end_ns),
                text_param(outcome.as_str()),
                opt_text_param(diag),
                SqliteValue::Integer(op),
            ])
            .map_err(|e| sql_err("op finish", &e))?;
        if affected == 1 {
            return Ok(());
        }
        match self.op(op)? {
            Some(_) => Err(LedgerError::DoubleFinish { op }),
            None => Err(LedgerError::NotFound {
                what: format!("op {op}"),
            }),
        }
    }

    fn op_row_is_bounded(&self, id: i64) -> Result<bool, LedgerError> {
        let rows = self
            .conn
            .query_with_params(
                &format!(
                    "SELECT CASE WHEN {} THEN 1 ELSE 0 END FROM ops WHERE id = ?1",
                    op_storage_predicate()
                ),
                &[SqliteValue::Integer(id)],
            )
            .map_err(|e| sql_err("op bounded-storage preflight", &e))?;
        let Some(row) = rows.first() else {
            return Ok(false);
        };
        if row_i64(row, 0, "op bounded-storage preflight")? != 1 {
            return Err(LedgerError::OpCorrupt {
                op: id,
                detail: "a variable-size field has the wrong storage type, exceeds its byte \
                         bound, contains invalid JSON, or the outcome envelope is malformed"
                    .to_string(),
            });
        }
        Ok(true)
    }

    fn op_execution_context_inner(
        &self,
        id: i64,
    ) -> Result<Option<OpExecutionContext>, LedgerError> {
        if !self.op_row_is_bounded(id)? {
            return Ok(None);
        }
        let rows = self
            .conn
            .query_with_params(
                &format!(
                    "SELECT branch, exec_mode, \
                            EXISTS(SELECT 1 FROM branches WHERE id = ops.branch) \
                     FROM ops WHERE id = ?1 AND {}",
                    op_storage_predicate()
                ),
                &[SqliteValue::Integer(id)],
            )
            .map_err(|e| sql_err("bounded op execution context", &e))?;
        let row = rows.first().ok_or_else(|| LedgerError::OpCorrupt {
            op: id,
            detail: "operation context disappeared or changed after bounded preflight".to_string(),
        })?;
        let branch = match row.get(0) {
            Some(SqliteValue::Integer(branch)) => *branch,
            _ => {
                return Err(LedgerError::OpCorrupt {
                    op: id,
                    detail: "branch has non-integer storage after bounded preflight".to_string(),
                });
            }
        };
        let exec_mode = match row.get(1) {
            Some(SqliteValue::Text(mode)) => {
                ExecMode::parse(mode).ok_or_else(|| LedgerError::OpCorrupt {
                    op: id,
                    detail: "execution mode is outside the deterministic/fast domain".to_string(),
                })?
            }
            _ => {
                return Err(LedgerError::OpCorrupt {
                    op: id,
                    detail: "execution mode has non-text storage after bounded preflight"
                        .to_string(),
                });
            }
        };
        if !matches!(row.get(2), Some(SqliteValue::Integer(1))) {
            return Err(LedgerError::OpCorrupt {
                op: id,
                detail: format!("branch {branch} does not exist in the branch registry"),
            });
        }
        Ok(Some(OpExecutionContext { branch, exec_mode }))
    }

    /// Fetch the fixed-size branch/mode execution context for one op.
    ///
    /// The complete bounded op envelope is checked before this query returns,
    /// but its variable-size fields are never materialized.
    ///
    /// # Errors
    /// Engine errors or [`LedgerError::OpCorrupt`] for a malformed operation
    /// or missing branch; an unknown op is `Ok(None)`.
    pub fn op_execution_context(&self, id: i64) -> Result<Option<OpExecutionContext>, LedgerError> {
        self.note_read_query();
        self.op_execution_context_inner(id)
    }

    pub(crate) fn bounded_op_exec_mode(&self, id: i64) -> Result<String, LedgerError> {
        self.op_execution_context_inner(id)?
            .map(|context| context.exec_mode.as_str().to_string())
            .ok_or_else(|| LedgerError::NotFound {
                what: format!("op {id}"),
            })
    }

    /// Fetch one op row, if present. Every variable-size field is checked by
    /// a metadata-only SQL preflight before the guarded payload query can
    /// materialize it.
    ///
    /// # Errors
    /// Engine errors or [`LedgerError::OpCorrupt`] when a stored envelope
    /// violates the bounded read contract; absence is `Ok(None)`.
    pub fn op(&self, id: i64) -> Result<Option<OpRow>, LedgerError> {
        self.note_read_query();
        if !self.op_row_is_bounded(id)? {
            return Ok(None);
        }
        let rows = self
            .conn
            .query_with_params(
                &format!(
                    "SELECT id, session, ir, seed, versions, budget, capability, t_start, t_end, \
                     outcome, diag FROM ops WHERE id = ?1 AND {}",
                    op_storage_predicate()
                ),
                &[SqliteValue::Integer(id)],
            )
            .map_err(|e| sql_err("op fetch", &e))?;
        let row = rows.first().ok_or_else(|| LedgerError::OpCorrupt {
            op: id,
            detail: "operation row disappeared or exceeded a guarded field bound after preflight"
                .to_string(),
        })?;
        let text_at = |idx: usize| -> Result<String, LedgerError> {
            match row.get(idx) {
                Some(SqliteValue::Text(t)) => Ok(t.as_str().to_string()),
                other => Err(LedgerError::Sql {
                    context: "op fetch".to_string(),
                    detail: format!("column {idx}: expected TEXT, got {other:?}"),
                }),
            }
        };
        let opt_text_at = |idx: usize| -> Option<String> {
            match row.get(idx) {
                Some(SqliteValue::Text(t)) => Some(t.as_str().to_string()),
                _ => None,
            }
        };
        let session = match row.get(1) {
            Some(SqliteValue::Blob(b)) => Some(b.to_vec()),
            _ => None,
        };
        let seed = match row.get(3) {
            Some(SqliteValue::Blob(b)) => b.to_vec(),
            other => {
                return Err(LedgerError::Sql {
                    context: "op fetch".to_string(),
                    detail: format!("seed: expected BLOB, got {other:?}"),
                });
            }
        };
        let t_end = match row.get(8) {
            Some(SqliteValue::Integer(v)) => Some(*v),
            _ => None,
        };
        Ok(Some(OpRow {
            id: row_i64(row, 0, "op.id")?,
            session,
            ir: text_at(2)?,
            seed,
            versions: text_at(4)?,
            budget: text_at(5)?,
            capability: text_at(6)?,
            t_start: row_i64(row, 7, "op.t_start")?,
            t_end,
            outcome: opt_text_at(9),
            diag: opt_text_at(10),
        }))
    }

    fn zero_row_link_refusal(&self, op: i64) -> Result<(), LedgerError> {
        let Some(_) = self.op_execution_context_inner(op)? else {
            return Err(LedgerError::NotFound {
                what: format!("op {op}"),
            });
        };
        let rows = self
            .conn
            .query_with_params(
                "SELECT CASE WHEN t_end IS NULL AND outcome IS NULL THEN 0 ELSE 1 END \
                 FROM ops WHERE id = ?1 LIMIT 1",
                &[SqliteValue::Integer(op)],
            )
            .map_err(|error| sql_err("classify refused edge insert", &error))?;
        let row = rows.first().ok_or_else(|| LedgerError::NotFound {
            what: format!("op {op}"),
        })?;
        match row_i64(row, 0, "classify refused edge insert")? {
            1 => Err(LedgerError::OpLineageSealed { op }),
            0 => Err(LedgerError::OpCorrupt {
                op,
                detail: "the atomic edge insert selected zero rows even though the operation \
                         remained bounded, branch-valid, and unfinished"
                    .to_string(),
            }),
            _ => Err(LedgerError::OpCorrupt {
                op,
                detail: "the fixed-size finish-state classifier returned a non-boolean value"
                    .to_string(),
            }),
        }
    }

    /// Link an unfinished op to an artifact in the lineage DAG. The guarded
    /// `INSERT ... SELECT` and [`Ledger::finish_op`] are serialized by the
    /// database: the link either commits before the finish, or changes zero
    /// rows and returns [`LedgerError::OpLineageSealed`] after it. Foreign keys
    /// still require the artifact to exist.
    ///
    /// # Errors
    /// [`LedgerError::NotFound`] for an unknown op;
    /// [`LedgerError::OpLineageSealed`] for a finished op;
    /// [`LedgerError::OpCorrupt`] for a malformed or branch-orphaned op;
    /// [`LedgerError::Invalid`] for a missing artifact or a conflicting
    /// explicit lineage seal.
    pub fn link(&self, op: i64, artifact: &ContentHash, role: EdgeRole) -> Result<(), LedgerError> {
        let bounded_op = op_storage_predicate();
        let insert = self
            .conn
            .prepare(&format!(
                "INSERT INTO edges(op, artifact, role) \
                 SELECT ops.id, ?2, ?3 FROM ops \
                 WHERE ops.id = ?1 AND ops.t_end IS NULL AND ops.outcome IS NULL \
                   AND ({bounded_op}) \
                   AND EXISTS(SELECT 1 FROM branches WHERE branches.id = ops.branch)"
            ))
            .map_err(|e| sql_err("edge insert prepare", &e))?
            .execute_with_params(&[
                SqliteValue::Integer(op),
                blob_param(artifact.as_bytes()),
                text_param(role.as_str()),
            ]);
        match insert {
            Ok(1) => Ok(()),
            Ok(0) => self.zero_row_link_refusal(op),
            Ok(changed) => Err(LedgerError::OpCorrupt {
                op,
                detail: format!(
                    "one edge insert changed {changed} rows; exactly one row was required"
                ),
            }),
            Err(FrankenError::ForeignKeyViolation) => {
                match self.artifact_content_identity(artifact) {
                    Err(error) => Err(error),
                    Ok(None) => Err(LedgerError::Invalid {
                        field: "edge".to_string(),
                        problem: format!(
                            "artifact {} does not exist; record it before linking to op {op}",
                            artifact.to_hex()
                        ),
                    }),
                    Ok(Some(_)) => Err(LedgerError::OpCorrupt {
                        op,
                        detail: "edge insertion hit a foreign-key refusal even though the linked artifact and typed content identity both reverified"
                            .to_string(),
                    }),
                }
            }
            Err(error)
                if error
                    .to_string()
                    .contains("sealed artifact rejects a different output producer") =>
            {
                Err(LedgerError::Invalid {
                    field: "edge".to_string(),
                    problem: format!(
                        "artifact {} has an immutable exclusive output-producer seal",
                        artifact.to_hex()
                    ),
                })
            }
            Err(error)
                if error
                    .to_string()
                    .contains("sealed operation artifact-edge set is immutable") =>
            {
                Err(LedgerError::Invalid {
                    field: "edge".to_string(),
                    problem: format!(
                        "operation {op} has an immutable artifact-edge-set seal; no edge may be added"
                    ),
                })
            }
            Err(e) => Err(sql_err("edge insert", &e)),
        }
    }

    /// Whether the lineage DAG contains this exact role-qualified edge.
    ///
    /// This is the verifier-side companion to [`Ledger::link`]: callers can
    /// prove that a content-addressed artifact was an input or output of the
    /// claimed operation without scanning or reconstructing the whole DAG.
    ///
    /// # Errors
    /// Engine errors only. Missing operations, artifacts, or edges return
    /// `Ok(false)`.
    pub fn edge_exists(
        &self,
        op: i64,
        artifact: &ContentHash,
        role: EdgeRole,
    ) -> Result<bool, LedgerError> {
        self.note_read_query();
        let rows = self
            .conn
            .query_with_params(
                "SELECT 1 FROM edges WHERE op = ?1 AND artifact = ?2 AND role = ?3 LIMIT 1",
                &[
                    SqliteValue::Integer(op),
                    blob_param(artifact.as_bytes()),
                    text_param(role.as_str()),
                ],
            )
            .map_err(|error| sql_err("edge existence query", &error))?;
        Ok(!rows.is_empty())
    }

    fn artifact_output_seal_inner(
        &self,
        artifact: &ContentHash,
    ) -> Result<Option<i64>, LedgerError> {
        let rows = self
            .conn
            .query_with_params(
                "SELECT CASE WHEN typeof(seal.op) = 'integer' THEN seal.op ELSE NULL END, \
                        CASE WHEN typeof(seal.role) = 'text' AND seal.role = 'out' \
                             THEN 1 ELSE 0 END, \
                        EXISTS( \
                            SELECT 1 FROM edges INDEXED BY idx_edges_artifact_role_op \
                            WHERE artifact = seal.artifact AND role = 'out' AND op = seal.op \
                            LIMIT 1 \
                        ), \
                        NOT EXISTS( \
                            SELECT 1 FROM edges INDEXED BY idx_edges_artifact_role_op \
                            WHERE artifact = seal.artifact AND role = 'out' AND op != seal.op \
                            LIMIT 1 \
                        ) \
                 FROM artifact_output_seals AS seal WHERE seal.artifact = ?1 LIMIT 1",
                &[blob_param(artifact.as_bytes())],
            )
            .map_err(|error| sql_err("artifact output seal query", &error))?;
        let Some(row) = rows.first() else {
            return Ok(None);
        };
        let op = match row.get(0) {
            Some(SqliteValue::Integer(op)) => *op,
            _ => {
                return Err(LedgerError::Corrupt {
                    hash_hex: artifact.to_hex(),
                    detail: "artifact output seal has a non-integer operation identity".to_string(),
                });
            }
        };
        if !matches!(row.get(1), Some(SqliteValue::Integer(1)))
            || !matches!(row.get(2), Some(SqliteValue::Integer(1)))
            || !matches!(row.get(3), Some(SqliteValue::Integer(1)))
        {
            return Err(LedgerError::Corrupt {
                hash_hex: artifact.to_hex(),
                detail: "artifact output seal has a malformed role, missing exact output edge, or competing producer"
                    .to_string(),
            });
        }
        Ok(Some(op))
    }

    /// Return the immutable exclusive output-producer seal for one artifact.
    ///
    /// A seal is absent unless a consumer explicitly requested single-producer
    /// provenance. Once present, schema-attested triggers reject every output
    /// edge from a different operation.
    ///
    /// # Errors
    /// Engine errors or [`LedgerError::Corrupt`] for malformed stored state.
    pub fn artifact_output_seal(&self, artifact: &ContentHash) -> Result<Option<i64>, LedgerError> {
        self.note_read_query();
        self.artifact_output_seal_inner(artifact)
    }

    /// Atomically seal one artifact to its existing sole output producer.
    ///
    /// The insert trigger verifies that `(op, artifact, out)` exists and that
    /// no different output producer exists in the same SQLite statement. The
    /// seal is idempotent for the same `op`, immutable, and may be created
    /// inside a caller-owned transaction so op/edge/seal materialization can
    /// commit together.
    ///
    /// # Errors
    /// [`LedgerError::Invalid`] when the requested producer conflicts with an
    /// existing seal or the artifact does not currently have exactly that sole
    /// output producer; engine errors otherwise.
    pub fn seal_artifact_output(&self, artifact: &ContentHash, op: i64) -> Result<(), LedgerError> {
        match self.artifact_output_seal_inner(artifact)? {
            Some(stored) if stored == op => return Ok(()),
            Some(stored) => {
                return Err(LedgerError::Invalid {
                    field: "artifact_output_seal".to_string(),
                    problem: format!(
                        "artifact {} is already sealed to operation {stored}, not {op}",
                        artifact.to_hex()
                    ),
                });
            }
            None => {}
        }

        // The exact-producer predicate must live in trg_artifact_output_seals_exact_producer,
        // not in a guarded `INSERT ... SELECT ?1 WHERE EXISTS(... ?1 ...)`: fsqlite never
        // binds parameters referenced inside WHERE-clause subqueries of a FROM-less SELECT,
        // so that guard silently inserts zero rows even for a valid sole producer
        // (bead frankensim-lnbzs). The trigger's NEW.-based subqueries evaluate correctly
        // and abort in the same statement, preserving atomicity.
        let insert = self
            .conn
            .prepare("INSERT INTO artifact_output_seals(artifact, op, role) VALUES (?1, ?2, 'out')")
            .map_err(|error| sql_err("artifact output seal prepare", &error))?
            .execute_with_params(&[blob_param(artifact.as_bytes()), SqliteValue::Integer(op)]);
        match insert {
            Ok(1) => Ok(()),
            Ok(0) => Err(LedgerError::Invalid {
                field: "artifact_output_seal".to_string(),
                problem: format!(
                    "artifact {} does not have operation {op} as its sole output producer",
                    artifact.to_hex()
                ),
            }),
            Ok(affected) => Err(LedgerError::Sql {
                context: "artifact output seal insert".to_string(),
                detail: format!("expected one inserted row, observed {affected}"),
            }),
            Err(error) => match self.artifact_output_seal_inner(artifact)? {
                Some(stored) if stored == op => Ok(()),
                Some(stored) => Err(LedgerError::Invalid {
                    field: "artifact_output_seal".to_string(),
                    problem: format!(
                        "artifact {} raced with an immutable seal for operation {stored}, not {op}",
                        artifact.to_hex()
                    ),
                }),
                None if error
                    .to_string()
                    .contains("artifact output seal requires one exact producer") =>
                {
                    Err(LedgerError::Invalid {
                        field: "artifact_output_seal".to_string(),
                        problem: format!(
                            "artifact {} does not have operation {op} as its sole output producer",
                            artifact.to_hex()
                        ),
                    })
                }
                None => Err(sql_err("artifact output seal insert", &error)),
            },
        }
    }

    fn op_artifact_edge_seal_inner(&self, op: i64) -> Result<Option<usize>, LedgerError> {
        let rows = self
            .conn
            .query_with_params(
                &format!(
                    "SELECT \
                         CASE WHEN typeof(seal.edge_count) = 'integer' AND \
                                        seal.edge_count BETWEEN 0 AND {MAX_LINEAGE_QUERY_ROWS} \
                              THEN seal.edge_count ELSE NULL END, \
                         EXISTS(SELECT 1 FROM ops WHERE id = seal.op LIMIT 1) \
                     FROM op_artifact_edge_seals AS seal WHERE seal.op = ?1 LIMIT 1"
                ),
                &[SqliteValue::Integer(op)],
            )
            .map_err(|error| sql_err("op artifact-edge seal query", &error))?;
        let Some(row) = rows.first() else {
            return Ok(None);
        };
        let edge_count = match row.get(0) {
            Some(SqliteValue::Integer(count)) => {
                usize::try_from(*count).map_err(|_| LedgerError::OpCorrupt {
                    op,
                    detail: "artifact-edge seal count is negative or unrepresentable".to_string(),
                })?
            }
            _ => {
                return Err(LedgerError::OpCorrupt {
                    op,
                    detail: "artifact-edge seal has a non-integer or out-of-range count"
                        .to_string(),
                });
            }
        };
        if !matches!(row.get(1), Some(SqliteValue::Integer(1))) {
            return Err(LedgerError::OpCorrupt {
                op,
                detail: "artifact-edge seal references a missing operation".to_string(),
            });
        }
        let probe_limit = edge_count + 1;
        let actual = self
            .conn
            .query_with_params(
                &format!(
                    "SELECT 1 FROM edges INDEXED BY idx_edges_op_role_artifact \
                     WHERE op = ?1 LIMIT {probe_limit}"
                ),
                &[SqliteValue::Integer(op)],
            )
            .map_err(|error| sql_err("op artifact-edge seal validation", &error))?
            .len();
        if actual != edge_count {
            return Err(LedgerError::OpCorrupt {
                op,
                detail: format!(
                    "artifact-edge seal records {edge_count} edges but the bounded validation observed {actual}"
                ),
            });
        }
        Ok(Some(edge_count))
    }

    /// Return the immutable artifact-edge-set seal for one operation.
    ///
    /// The stored count and operation parent are fixed-size validated, then a
    /// covering-index probe reads at most `count + 1` rows to prove the sealed
    /// set has not been bypass-mutated.
    ///
    /// # Errors
    /// Engine errors or [`LedgerError::OpCorrupt`] for malformed or
    /// count-inconsistent stored state.
    pub fn op_artifact_edge_seal(&self, op: i64) -> Result<Option<usize>, LedgerError> {
        self.note_read_query();
        self.op_artifact_edge_seal_inner(op)
    }

    /// Atomically freeze one operation's complete current artifact-edge set.
    ///
    /// `expected_count` is capped at [`MAX_LINEAGE_QUERY_ROWS`]. The insert
    /// reads at most `expected_count + 1` covering-index rows and succeeds only
    /// on exact cardinality; schema-attested triggers then reject every edge
    /// insert, update, or delete involving the sealed operation. Exact same-
    /// count retry is idempotent.
    ///
    /// # Errors
    /// [`LedgerError::Invalid`] for an excessive cap, missing op, cardinality
    /// mismatch, or conflicting existing seal; engine errors otherwise.
    pub fn seal_op_artifact_edges(
        &self,
        op: i64,
        expected_count: usize,
    ) -> Result<(), LedgerError> {
        if expected_count > MAX_LINEAGE_QUERY_ROWS {
            return Err(LedgerError::Invalid {
                field: "op_artifact_edge_seal".to_string(),
                problem: format!(
                    "edge count {expected_count} exceeds the public maximum {MAX_LINEAGE_QUERY_ROWS}"
                ),
            });
        }
        match self.op_artifact_edge_seal_inner(op)? {
            Some(stored) if stored == expected_count => return Ok(()),
            Some(stored) => {
                return Err(LedgerError::Invalid {
                    field: "op_artifact_edge_seal".to_string(),
                    problem: format!(
                        "operation {op} is already sealed at {stored} edges, not {expected_count}"
                    ),
                });
            }
            None => {}
        }

        let expected_i64 = i64::try_from(expected_count).expect("public lineage cap fits i64");
        // Same constraint as seal_artifact_output: the op-exists + exact-bounded-count
        // predicate must live in trg_op_artifact_edge_seals_exact_count because fsqlite
        // never binds parameters inside WHERE-clause subqueries of a FROM-less
        // `INSERT ... SELECT` guard (bead frankensim-lnbzs). The trigger bounds its
        // count probe at LIMIT 1025 (MAX_LINEAGE_QUERY_ROWS + 1) and aborts in the
        // same statement.
        let insert = self
            .conn
            .prepare("INSERT INTO op_artifact_edge_seals(op, edge_count) VALUES (?1, ?2)")
            .map_err(|error| sql_err("op artifact-edge seal prepare", &error))?
            .execute_with_params(&[SqliteValue::Integer(op), SqliteValue::Integer(expected_i64)]);
        match insert {
            Ok(1) => Ok(()),
            Ok(0) => Err(LedgerError::Invalid {
                field: "op_artifact_edge_seal".to_string(),
                problem: format!(
                    "operation {op} does not exist or does not have exactly {expected_count} artifact edges"
                ),
            }),
            Ok(affected) => Err(LedgerError::Sql {
                context: "op artifact-edge seal insert".to_string(),
                detail: format!("expected one inserted row, observed {affected}"),
            }),
            Err(error) => match self.op_artifact_edge_seal_inner(op)? {
                Some(stored) if stored == expected_count => Ok(()),
                Some(stored) => Err(LedgerError::Invalid {
                    field: "op_artifact_edge_seal".to_string(),
                    problem: format!(
                        "operation {op} raced with an immutable {stored}-edge seal, not {expected_count}"
                    ),
                }),
                None if error
                    .to_string()
                    .contains("op artifact-edge seal requires the exact bounded edge count") =>
                {
                    Err(LedgerError::Invalid {
                        field: "op_artifact_edge_seal".to_string(),
                        problem: format!(
                            "operation {op} does not have exactly {expected_count} artifact edges"
                        ),
                    })
                }
                None => Err(sql_err("op artifact-edge seal insert", &error)),
            },
        }
    }

    /// Return output-producer op ids for one artifact under an explicit row
    /// cap. The query reads at most `cap + 1` fixed-size rows; callers must
    /// reject or otherwise account for [`BoundedProducerOps::truncated`].
    ///
    /// A zero cap is valid and acts as a bounded existence probe. Unknown
    /// artifacts return an empty, non-truncated result.
    ///
    /// # Errors
    /// [`LedgerError::Invalid`] when `cap` exceeds
    /// [`MAX_LINEAGE_QUERY_ROWS`], [`LedgerError::Corrupt`] for malformed
    /// stored producer identities, or an engine error.
    pub fn artifact_producer_ops_bounded(
        &self,
        artifact: &ContentHash,
        cap: usize,
    ) -> Result<BoundedProducerOps, LedgerError> {
        let limit = bounded_lineage_query_limit(cap)?;
        self.note_read_query();
        let rows = self
            .conn
            .query_with_params(
                &format!(
                    "SELECT CASE WHEN typeof(op) = 'integer' THEN op ELSE NULL END \
                     FROM edges INDEXED BY idx_edges_artifact_role_op \
                     WHERE artifact = ?1 AND role = 'out' \
                     ORDER BY op LIMIT {limit}"
                ),
                &[blob_param(artifact.as_bytes())],
            )
            .map_err(|error| sql_err("bounded artifact producer query", &error))?;
        let truncated = rows.len() > cap;
        let mut op_ids = Vec::with_capacity(rows.len().min(cap));
        for row in rows.iter().take(cap) {
            match row.get(0) {
                Some(SqliteValue::Integer(op)) => op_ids.push(*op),
                _ => {
                    return Err(LedgerError::Corrupt {
                        hash_hex: artifact.to_hex(),
                        detail: "an output-producer edge has a non-integer operation identity"
                            .to_string(),
                    });
                }
            }
        }
        Ok(BoundedProducerOps { op_ids, truncated })
    }

    /// Return role-qualified artifact edges for one op under an explicit row
    /// cap. The query reads at most `cap + 1` fixed-size rows and orders them
    /// deterministically by role then content identity.
    ///
    /// A zero cap is valid and acts as a bounded existence probe. Unknown
    /// operations return an empty, non-truncated result.
    ///
    /// # Errors
    /// [`LedgerError::Invalid`] when `cap` exceeds
    /// [`MAX_LINEAGE_QUERY_ROWS`], [`LedgerError::OpCorrupt`] for malformed
    /// stored edges, or an engine error.
    pub fn op_artifact_edges_bounded(
        &self,
        op: i64,
        cap: usize,
    ) -> Result<BoundedOpArtifactEdges, LedgerError> {
        let limit = bounded_lineage_query_limit(cap)?;
        self.note_read_query();
        let rows = self
            .conn
            .query_with_params(
                &format!(
                    "SELECT \
                         CASE WHEN typeof(role) = 'text' AND role IN ('in','out') \
                              THEN role ELSE NULL END, \
                         CASE WHEN typeof(artifact) = 'blob' AND length(artifact) = 32 \
                              THEN artifact ELSE NULL END \
                     FROM edges INDEXED BY idx_edges_op_role_artifact WHERE op = ?1 \
                     ORDER BY role, artifact LIMIT {limit}"
                ),
                &[SqliteValue::Integer(op)],
            )
            .map_err(|error| sql_err("bounded op artifact-edge query", &error))?;
        let truncated = rows.len() > cap;
        let mut edges = Vec::with_capacity(rows.len().min(cap));
        for row in rows.iter().take(cap) {
            let role = match row.get(0) {
                Some(SqliteValue::Text(role)) => {
                    EdgeRole::parse(role).ok_or_else(|| LedgerError::OpCorrupt {
                        op,
                        detail: "an artifact edge has a role outside the in/out domain".to_string(),
                    })?
                }
                _ => {
                    return Err(LedgerError::OpCorrupt {
                        op,
                        detail: "an artifact edge role has non-text storage".to_string(),
                    });
                }
            };
            let artifact = match row.get(1) {
                Some(SqliteValue::Blob(bytes)) => {
                    ContentHash::from_slice(bytes).ok_or_else(|| LedgerError::OpCorrupt {
                        op,
                        detail: "an artifact edge has a malformed content identity".to_string(),
                    })?
                }
                _ => {
                    return Err(LedgerError::OpCorrupt {
                        op,
                        detail: "an artifact edge identity has non-blob storage".to_string(),
                    });
                }
            };
            edges.push(OpArtifactEdge { role, artifact });
        }
        Ok(BoundedOpArtifactEdges { edges, truncated })
    }

    // -- metrics, events, tune ---------------------------------------------

    /// Append one metric sample. `value` must be finite (REAL NOT NULL).
    ///
    /// # Errors
    /// [`LedgerError::Invalid`] for non-finite values; engine errors
    /// otherwise (including PK conflicts on duplicate `(op, t, name)`).
    pub fn record_metric(
        &self,
        op: i64,
        t: i64,
        name: &str,
        value: f64,
    ) -> Result<(), LedgerError> {
        if !value.is_finite() {
            return Err(LedgerError::Invalid {
                field: "value".to_string(),
                problem: format!("{value} is not finite; metrics are REAL NOT NULL"),
            });
        }
        self.conn
            .prepare("INSERT INTO metrics(op, t, name, value) VALUES (?1, ?2, ?3, ?4)")
            .map_err(|e| sql_err("metric insert prepare", &e))?
            .execute_with_params(&[
                SqliteValue::Integer(op),
                SqliteValue::Integer(t),
                text_param(name),
                SqliteValue::Float(value),
            ])
            .map_err(|e| sql_err("metric insert", &e))?;
        Ok(())
    }

    /// Append one event-stream row; returns its rowid.
    ///
    /// # Errors
    /// [`LedgerError::Invalid`] for malformed payload JSON or a reserved
    /// internal event kind; engine errors.
    pub fn append_event(&self, event: &EventRow<'_>) -> Result<i64, LedgerError> {
        if event.kind == VCS_IDENTITY_EVENT_KIND {
            return Err(LedgerError::Invalid {
                field: "kind".to_string(),
                problem: format!(
                    "event kind {VCS_IDENTITY_EVENT_KIND:?} is reserved for ledger identity"
                ),
            });
        }
        self.append_event_unchecked(event)
    }

    pub(crate) fn append_vcs_identity_event(
        &self,
        event: &EventRow<'_>,
    ) -> Result<i64, LedgerError> {
        if event.kind != VCS_IDENTITY_EVENT_KIND {
            return Err(LedgerError::Invalid {
                field: "kind".to_string(),
                problem: "the VCS identity insertion path accepts only its reserved event kind"
                    .to_string(),
            });
        }
        self.append_event_unchecked(event)
    }

    fn append_event_unchecked(&self, event: &EventRow<'_>) -> Result<i64, LedgerError> {
        if let Some(p) = event.payload {
            self.require_json("payload", p, false)?;
        }
        self.conn
            .prepare("INSERT INTO events(session, t, kind, payload) VALUES (?1, ?2, ?3, ?4)")
            .map_err(|e| sql_err("event insert prepare", &e))?
            .execute_with_params(&[
                opt_blob_param(event.session),
                SqliteValue::Integer(event.t),
                text_param(event.kind),
                opt_text_param(event.payload),
            ])
            .map_err(|e| sql_err("event insert", &e))?;
        Ok(self.conn.last_insert_rowid())
    }

    /// Append a batch of events in one transaction (the append-heavy write
    /// path of plan §11.2).
    ///
    /// # Errors
    /// On any failure the whole batch rolls back (when this call owns the
    /// transaction).
    pub fn append_events(&self, batch: &[EventRow<'_>]) -> Result<(), LedgerError> {
        let owns_txn = !self.conn.in_transaction();
        if owns_txn {
            self.begin()?;
        }
        for event in batch {
            if let Err(e) = self.append_event(event) {
                if owns_txn {
                    let _ = self.rollback();
                }
                return Err(e);
            }
        }
        if owns_txn && let Err(e) = self.commit() {
            let _ = self.rollback();
            return Err(e);
        }
        Ok(())
    }

    /// Number of stored events (stress/verification helper).
    ///
    /// # Errors
    /// Engine errors.
    pub fn table_count(&self, table: &str) -> Result<u64, LedgerError> {
        if !ALL_TABLES.contains(&table) {
            return Err(LedgerError::Invalid {
                field: "table".to_string(),
                problem: format!("unknown table {table:?}; see fs_ledger::ALL_TABLES"),
            });
        }
        let row = self
            .conn
            .query_row(&format!("SELECT COUNT(*) FROM {table}"))
            .map_err(|e| sql_err("table_count", &e))?;
        row_u64(&row, 0, "table_count")
    }

    /// Upsert one autotuner cache row (`kernel` × `shape_class` × machine
    /// fingerprint).
    ///
    /// The cache-row upsert remains one `INSERT .. ON CONFLICT .. DO UPDATE`;
    /// schema-v19 binaries publish its typed content sidecar in the same
    /// transaction. (Between 2026-07-11 and the upstream fix
    /// 3c388122 this was routed UPDATE-then-INSERT around the fsqlite
    /// upsert-after-leaf-split corruption, bead u8og /
    /// Dicklesworthstone/frankensqlite#123; the raw repro from that issue
    /// was re-verified clean on both ISAs before restoring this form.)
    ///
    /// # Errors
    /// [`LedgerError::Invalid`] for a non-canonical or oversized key/blob/JSON;
    /// engine errors otherwise.
    pub fn tune_put(
        &self,
        kernel: &str,
        shape_class: &str,
        machine: &[u8],
        params: &str,
        measured: &str,
    ) -> Result<(), LedgerError> {
        self.require_tune_row(kernel, shape_class, machine, params, measured)?;
        let owns_txn = !self.conn.in_transaction();
        if owns_txn {
            self.begin()?;
        }
        let result = (|| {
            let changed = self
                .conn
                .prepare(
                    "INSERT INTO tune(kernel, shape_class, machine, params, measured) \
                     VALUES (?1, ?2, ?3, ?4, ?5) \
                     ON CONFLICT(kernel, shape_class, machine) \
                     DO UPDATE SET params = excluded.params, measured = excluded.measured",
                )
                .map_err(|e| sql_err("tune upsert prepare", &e))?
                .execute_with_params(&[
                    text_param(kernel),
                    text_param(shape_class),
                    blob_param(machine),
                    text_param(params),
                    text_param(measured),
                ])
                .map_err(|e| sql_err("tune upsert", &e))?;
            if changed != 1 {
                return Err(tune_corrupt(
                    kernel,
                    format!("one tune upsert changed {changed} rows"),
                ));
            }
            self.upsert_tune_content_identity(kernel, shape_class, machine, params, measured)
        })();
        match (&result, owns_txn) {
            (Ok(()), true) => {
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
        result
    }

    /// Insert one autotuner row only when its exact storage key is absent.
    /// Existing rows are never modified. Callers that require idempotent exact
    /// identity should fetch and compare after this call.
    ///
    /// # Errors
    /// [`LedgerError::Invalid`] for a non-canonical or oversized key/blob/JSON;
    /// engine errors otherwise.
    pub fn tune_put_if_absent(
        &self,
        kernel: &str,
        shape_class: &str,
        machine: &[u8],
        params: &str,
        measured: &str,
    ) -> Result<(), LedgerError> {
        self.require_tune_row(kernel, shape_class, machine, params, measured)?;
        let owns_txn = !self.conn.in_transaction();
        if owns_txn {
            self.begin()?;
        }
        let result = (|| {
            let changed = self
                .conn
                .prepare(
                    "INSERT INTO tune(kernel, shape_class, machine, params, measured) \
                     VALUES (?1, ?2, ?3, ?4, ?5) \
                     ON CONFLICT(kernel, shape_class, machine) DO NOTHING",
                )
                .map_err(|e| sql_err("tune insert-if-absent prepare", &e))?
                .execute_with_params(&[
                    text_param(kernel),
                    text_param(shape_class),
                    blob_param(machine),
                    text_param(params),
                    text_param(measured),
                ])
                .map_err(|e| sql_err("tune insert-if-absent", &e))?;
            match changed {
                1 => self.upsert_tune_content_identity(
                    kernel,
                    shape_class,
                    machine,
                    params,
                    measured,
                ),
                0 => self.verify_tune_content_identity(kernel, shape_class, machine),
                other => Err(tune_corrupt(
                    kernel,
                    format!("one tune insert-if-absent changed {other} rows"),
                )),
            }
        })();
        match (&result, owns_txn) {
            (Ok(()), true) => {
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
        result
    }

    /// Fetch one autotuner cache row, if present.
    ///
    /// # Errors
    /// [`LedgerError::Invalid`] for a non-canonical lookup key,
    /// [`LedgerError::TuneCorrupt`] for a stored row outside the bounded tune
    /// contract, or engine errors. Absence is `Ok(None)`.
    pub fn tune_get(
        &self,
        kernel: &str,
        shape_class: &str,
        machine: &[u8],
    ) -> Result<Option<TuneRow>, LedgerError> {
        self.note_read_query();
        self.tune_get_inner(kernel, shape_class, machine)
    }

    fn tune_get_inner(
        &self,
        kernel: &str,
        shape_class: &str,
        machine: &[u8],
    ) -> Result<Option<TuneRow>, LedgerError> {
        validate_tune_identity("kernel", kernel, MAX_TUNE_KERNEL_BYTES)?;
        validate_tune_identity("shape_class", shape_class, MAX_TUNE_SHAPE_CLASS_BYTES)?;
        validate_tune_machine(machine)?;

        let metadata_sql = format!(
            "SELECT typeof(params), length(CAST(params AS BLOB)), \
                    CASE WHEN typeof(params) = 'text' THEN \
                        CASE WHEN length(CAST(params AS BLOB)) BETWEEN 1 AND {MAX_TUNE_PARAMS_BYTES} \
                            THEN json_valid(params) ELSE 0 END \
                    ELSE 0 END, \
                    typeof(measured), length(CAST(measured AS BLOB)), \
                    CASE WHEN typeof(measured) = 'text' THEN \
                        CASE WHEN length(CAST(measured AS BLOB)) BETWEEN 1 AND {MAX_TUNE_MEASURED_BYTES} \
                            THEN json_valid(measured) ELSE 0 END \
                    ELSE 0 END \
             FROM tune \
             WHERE kernel = ?1 AND shape_class = ?2 AND machine = ?3"
        );
        let metadata = self
            .conn
            .query_with_params(
                &metadata_sql,
                &[
                    text_param(kernel),
                    text_param(shape_class),
                    blob_param(machine),
                ],
            )
            .map_err(|e| sql_err("tune fetch metadata", &e))?;
        let Some(metadata_row) = metadata.first() else {
            return Ok(None);
        };
        tune_column_len(
            metadata_row,
            TuneColumnSpec {
                type_idx: 0,
                len_idx: 1,
                field: "params",
                expected_type: "text",
                min_bytes: 1,
                max_bytes: MAX_TUNE_PARAMS_BYTES,
            },
            kernel,
        )?;
        require_stored_json(metadata_row, 2, "params", kernel)?;
        tune_column_len(
            metadata_row,
            TuneColumnSpec {
                type_idx: 3,
                len_idx: 4,
                field: "measured",
                expected_type: "text",
                min_bytes: 1,
                max_bytes: MAX_TUNE_MEASURED_BYTES,
            },
            kernel,
        )?;
        require_stored_json(metadata_row, 5, "measured", kernel)?;

        let guarded_sql = format!(
            "SELECT params, measured FROM tune \
             WHERE kernel = ?1 AND shape_class = ?2 AND machine = ?3 AND \
                   typeof(params) = 'text' AND \
                   length(CAST(params AS BLOB)) BETWEEN 1 AND {MAX_TUNE_PARAMS_BYTES} AND \
                   CASE WHEN typeof(params) = 'text' THEN \
                       CASE WHEN length(CAST(params AS BLOB)) BETWEEN 1 AND {MAX_TUNE_PARAMS_BYTES} \
                           THEN json_valid(params) ELSE 0 END \
                       ELSE 0 END = 1 AND \
                   typeof(measured) = 'text' AND \
                   length(CAST(measured AS BLOB)) BETWEEN 1 AND {MAX_TUNE_MEASURED_BYTES} AND \
                   CASE WHEN typeof(measured) = 'text' THEN \
                       CASE WHEN length(CAST(measured AS BLOB)) BETWEEN 1 AND {MAX_TUNE_MEASURED_BYTES} \
                           THEN json_valid(measured) ELSE 0 END \
                       ELSE 0 END = 1 LIMIT 1"
        );
        let guarded = self
            .conn
            .query_with_params(
                &guarded_sql,
                &[
                    text_param(kernel),
                    text_param(shape_class),
                    blob_param(machine),
                ],
            )
            .map_err(|e| sql_err("tune guarded fetch", &e))?;
        let Some(row) = guarded.first() else {
            return Err(LedgerError::Busy {
                context: "tune guarded fetch".to_string(),
                detail: "row changed after bounded metadata preflight".to_string(),
            });
        };
        Ok(Some(TuneRow {
            kernel: kernel.to_string(),
            shape_class: shape_class.to_string(),
            machine: machine.to_vec(),
            params: row_text(row, 0, "tune guarded fetch params")?,
            measured: row_text(row, 1, "tune guarded fetch measured")?,
        }))
    }

    /// All autotuner cache rows for one kernel, across shape classes and
    /// machine fingerprints (staleness scans: "a target that was never
    /// re-measured is a lie waiting to happen", plan §14.1).
    ///
    /// The result is ordered by `(shape_class, machine)`. The scan refuses
    /// before JSON/blob materialization when either the row cap or aggregate
    /// byte cap is exceeded.
    ///
    /// # Errors
    /// [`LedgerError::Invalid`] for a non-canonical kernel,
    /// [`LedgerError::TuneCorrupt`] for malformed stored rows,
    /// [`LedgerError::TuneReadLimit`] for a history outside the scan budget,
    /// or engine errors. An unknown kernel is an empty vec.
    pub fn tune_rows(&self, kernel: &str) -> Result<Vec<TuneRow>, LedgerError> {
        self.note_read_query();
        validate_tune_identity("kernel", kernel, MAX_TUNE_KERNEL_BYTES)?;
        let expected_rows = preflight_tune_scan(&self.conn, kernel)?;
        if expected_rows == 0 {
            return Ok(Vec::new());
        }
        guarded_tune_scan(&self.conn, kernel, expected_rows)
    }

    // -- Rev S extension tables ----------------------------------------------

    /// Upsert a named record in one of the Rev S extension tables.
    ///
    /// # Errors
    /// [`LedgerError::Invalid`] for empty name / malformed JSON; engine
    /// errors otherwise.
    pub fn put_extension(
        &self,
        table: ExtensionTable,
        name: &str,
        body_json: &str,
    ) -> Result<(), LedgerError> {
        if name.is_empty() {
            return Err(LedgerError::Invalid {
                field: "name".to_string(),
                problem: "empty; extension records are keyed by name".to_string(),
            });
        }
        self.require_json("body", body_json, false)?;
        self.conn
            .prepare(&format!(
                "INSERT INTO {}(name, body, created_at) VALUES (?1, ?2, ?3) \
                 ON CONFLICT(name) DO UPDATE SET body = excluded.body",
                table.table_name()
            ))
            .map_err(|e| sql_err("extension upsert prepare", &e))?
            .execute_with_params(&[
                text_param(name),
                text_param(body_json),
                SqliteValue::Integer(now_wall_ns()),
            ])
            .map_err(|e| sql_err("extension upsert", &e))?;
        Ok(())
    }

    /// Fetch a named extension record's JSON body, if present.
    ///
    /// # Errors
    /// Engine errors; absence is `Ok(None)`.
    pub fn get_extension(
        &self,
        table: ExtensionTable,
        name: &str,
    ) -> Result<Option<String>, LedgerError> {
        let rows = self
            .conn
            .query_with_params(
                &format!("SELECT body FROM {} WHERE name = ?1", table.table_name()),
                &[text_param(name)],
            )
            .map_err(|e| sql_err("extension fetch", &e))?;
        match rows.first().and_then(|r| r.get(0)) {
            Some(SqliteValue::Text(t)) => Ok(Some(t.as_str().to_string())),
            _ => Ok(None),
        }
    }

    // -- hygiene --------------------------------------------------------------

    /// Referential/shape hygiene scan. A crash-recovered ledger must lint
    /// clean: orphan edges are the acceptance criterion of the kill-storm
    /// battery.
    ///
    /// # Errors
    /// [`LedgerError::SchemaMismatch`] if an identity mutation guard is
    /// missing or changed; [`LedgerError::InstanceIdentityCorrupt`] if the
    /// singleton identity row is missing, malformed, or differs from this
    /// handle's cached authority; engine errors otherwise.
    #[allow(clippy::too_many_lines)] // One ordered, bounded cross-table hygiene report.
    pub fn lint(&self) -> Result<LintReport, LedgerError> {
        let _ = self.checked_instance_id()?;
        let count = |sql: &str, context: &str| -> Result<u64, LedgerError> {
            let row = self.conn.query_row(sql).map_err(|e| sql_err(context, &e))?;
            row_u64(&row, 0, context)
        };
        let malformed_storage_sql = format!(
            "SELECT COUNT(*) FROM artifacts a WHERE NOT \
             ((a.bytes IS NOT NULL AND a.chunk_count = 0) OR \
              (a.bytes IS NULL AND a.chunk_count > 0)) OR \
             typeof(a.kind) != 'text' OR \
             length(CAST(a.kind AS BLOB)) NOT BETWEEN 1 AND {MAX_ARTIFACT_KIND_BYTES} OR \
             (a.meta IS NOT NULL AND (typeof(a.meta) != 'text' OR \
              length(CAST(a.meta AS BLOB)) > {MAX_ARTIFACT_META_BYTES} OR \
              json_valid(a.meta) != 1)) OR \
             (a.bytes IS NOT NULL AND length(a.bytes) > {STORAGE_CHUNK_LEN}) OR \
             EXISTS (SELECT 1 FROM artifact_chunks c WHERE c.hash = a.hash \
                     AND (c.bytes IS NULL OR \
                          length(c.bytes) > {STORAGE_CHUNK_LEN}))"
        );
        let malformed_tune_sql = format!(
            "SELECT COUNT(*) FROM tune WHERE \
             typeof(kernel) != 'text' OR \
             length(CAST(kernel AS BLOB)) NOT BETWEEN 1 AND {MAX_TUNE_KERNEL_BYTES} OR \
             length(CAST(kernel AS BLOB)) != length(kernel) OR \
             kernel GLOB '*[^!-~]*' OR \
             typeof(shape_class) != 'text' OR \
             length(CAST(shape_class AS BLOB)) NOT BETWEEN 1 AND {MAX_TUNE_SHAPE_CLASS_BYTES} OR \
             length(CAST(shape_class AS BLOB)) != length(shape_class) OR \
             shape_class GLOB '*[^!-~]*' OR \
             typeof(machine) != 'blob' OR \
             length(machine) NOT BETWEEN 1 AND {MAX_TUNE_MACHINE_BYTES} OR \
             typeof(params) != 'text' OR \
             length(CAST(params AS BLOB)) NOT BETWEEN 1 AND {MAX_TUNE_PARAMS_BYTES} OR \
             CASE WHEN typeof(params) = 'text' THEN \
                 CASE WHEN length(CAST(params AS BLOB)) BETWEEN 1 AND {MAX_TUNE_PARAMS_BYTES} \
                     THEN json_valid(params) ELSE 0 END \
                 ELSE 0 END != 1 OR \
             typeof(measured) != 'text' OR \
             length(CAST(measured AS BLOB)) NOT BETWEEN 1 AND {MAX_TUNE_MEASURED_BYTES} OR \
             CASE WHEN typeof(measured) = 'text' THEN \
                 CASE WHEN length(CAST(measured AS BLOB)) BETWEEN 1 AND {MAX_TUNE_MEASURED_BYTES} \
                     THEN json_valid(measured) ELSE 0 END \
                 ELSE 0 END != 1"
        );
        let malformed_ops_sql = format!(
            "SELECT COUNT(*) FROM ops WHERE NOT ({})",
            op_storage_predicate()
        );
        Ok(LintReport {
            orphan_edge_ops: count(
                "SELECT COUNT(*) FROM edges e LEFT JOIN ops o ON e.op = o.id WHERE o.id IS NULL",
                "lint orphan_edge_ops",
            )?,
            orphan_edge_artifacts: count(
                "SELECT COUNT(*) FROM edges e LEFT JOIN artifacts a ON e.artifact = a.hash \
                 WHERE a.hash IS NULL",
                "lint orphan_edge_artifacts",
            )?,
            orphan_metric_ops: count(
                "SELECT COUNT(*) FROM metrics m LEFT JOIN ops o ON m.op = o.id \
                 WHERE o.id IS NULL",
                "lint orphan_metric_ops",
            )?,
            malformed_artifacts: count(&malformed_storage_sql, "lint malformed_artifacts")?,
            chunk_count_mismatches: count(
                "SELECT COUNT(*) FROM artifacts a WHERE a.chunk_count > 0 AND \
                 ((SELECT COUNT(*) FROM artifact_chunks c WHERE c.hash = a.hash) \
                    != a.chunk_count OR \
                  COALESCE((SELECT MIN(c.seq) FROM artifact_chunks c WHERE c.hash = a.hash), -1) \
                    != 0 OR \
                  COALESCE((SELECT MAX(c.seq) FROM artifact_chunks c WHERE c.hash = a.hash), -1) \
                    != a.chunk_count - 1)",
                "lint chunk_count_mismatches",
            )?,
            len_mismatches: count(
                "SELECT COUNT(*) FROM artifacts a WHERE \
                 (a.bytes IS NOT NULL AND length(a.bytes) != a.len) OR \
                 (a.chunk_count > 0 AND \
                  (SELECT COALESCE(SUM(length(c.bytes)), 0) FROM artifact_chunks c \
                   WHERE c.hash = a.hash) != a.len)",
                "lint len_mismatches",
            )?,
            orphan_chunks: count(
                "SELECT COUNT(*) FROM artifact_chunks c LEFT JOIN artifacts a \
                 ON c.hash = a.hash WHERE a.hash IS NULL",
                "lint orphan_chunks",
            )?,
            malformed_tune_rows: count(&malformed_tune_sql, "lint malformed_tune_rows")?,
            malformed_ops: count(&malformed_ops_sql, "lint malformed_ops")?,
            // AND/OR form rather than `(a IS NULL) != (b IS NULL)`: fsqlite
            // mis-associates postfix IS NULL against comparison operators
            // when re-parsing stored CHECK text (upstream bug; see bead).
            half_finished_ops: count(
                "SELECT COUNT(*) FROM ops WHERE \
                 (t_end IS NULL AND outcome IS NOT NULL) OR \
                 (t_end IS NOT NULL AND outcome IS NULL)",
                "lint half_finished_ops",
            )?,
            orphan_op_branches: count(
                "SELECT COUNT(*) FROM ops o LEFT JOIN branches b ON o.branch = b.id \
                 WHERE b.id IS NULL",
                "lint orphan_op_branches",
            )?,
            orphan_branch_parents: count(
                "SELECT COUNT(*) FROM branches c LEFT JOIN branches p ON c.parent = p.id \
                 WHERE c.parent IS NOT NULL AND p.id IS NULL",
                "lint orphan_branch_parents",
            )?,
            orphan_qty_crosswalk_sources: count(
                "SELECT COUNT(*) FROM qty_dimension_crosswalks q \
                 LEFT JOIN artifacts a ON q.old_hash = a.hash WHERE a.hash IS NULL",
                "lint orphan_qty_crosswalk_sources",
            )?,
            orphan_qty_crosswalk_targets: count(
                "SELECT COUNT(*) FROM qty_dimension_crosswalks q \
                 LEFT JOIN artifacts a ON q.new_hash = a.hash WHERE a.hash IS NULL",
                "lint orphan_qty_crosswalk_targets",
            )?,
            malformed_qty_crosswalks: count(
                "SELECT COUNT(*) FROM qty_dimension_crosswalks WHERE \
                 typeof(old_hash) != 'blob' OR length(old_hash) != 32 OR \
                 typeof(new_hash) != 'blob' OR length(new_hash) != 32 OR \
                 old_hash = new_hash OR \
                 typeof(source_version) != 'integer' OR source_version != 1 OR \
                 typeof(target_version) != 'integer' OR target_version != 2 OR \
                 typeof(rule) != 'text' OR rule != 'append-mole-zero' OR \
                 typeof(created_at) != 'integer'",
                "lint malformed_qty_crosswalks",
            )?,
            orphan_state_checkpoint_artifacts: count(
                "SELECT COUNT(*) FROM semantic_state_checkpoint_receipts c \
                 LEFT JOIN artifacts a ON c.runtime_state_artifact = a.hash \
                 WHERE a.hash IS NULL",
                "lint orphan_state_checkpoint_artifacts",
            )?,
            malformed_state_checkpoints: count(
                "SELECT COUNT(*) FROM semantic_state_checkpoint_receipts WHERE \
                 typeof(receipt_hash) != 'blob' OR length(receipt_hash) != 32 OR \
                 typeof(state_slot) != 'blob' OR length(state_slot) != 32 OR \
                 state_slot = X'0000000000000000000000000000000000000000000000000000000000000000' OR \
                 typeof(law_id) != 'blob' OR length(law_id) NOT BETWEEN 1 AND 256 OR \
                 typeof(law_version) != 'integer' OR \
                    law_version NOT BETWEEN 0 AND 4294967295 OR \
                 typeof(state_schema_version) != 'integer' OR \
                    state_schema_version NOT BETWEEN 0 AND 4294967295 OR \
                 typeof(runtime_state_artifact) != 'blob' OR \
                    length(runtime_state_artifact) != 32 OR \
                 typeof(canonical_parameters_hash) != 'blob' OR \
                    length(canonical_parameters_hash) != 32 OR \
                 canonical_parameters_hash = X'0000000000000000000000000000000000000000000000000000000000000000' OR \
                 typeof(contract_and_code_hash) != 'blob' OR \
                    length(contract_and_code_hash) != 32 OR \
                 contract_and_code_hash = X'0000000000000000000000000000000000000000000000000000000000000000' OR \
                 typeof(created_at) != 'integer'",
                "lint malformed_state_checkpoints",
            )?,
        })
    }
}

// ---------------------------------------------------------------------------
// Streaming artifact writer
// ---------------------------------------------------------------------------

static WRITER_NONCE: AtomicU64 = AtomicU64::new(0);

/// A provisional (non-content) chunk key for staging streamed chunks inside
/// the writer's transaction. Collision with a real BLAKE3 content hash would
/// require finding a preimage; treated as impossible (CONTRACT.md).
fn provisional_key() -> [u8; 32] {
    let mut h = Blake3::new();
    h.update(b"fs-ledger provisional chunk key v1");
    h.update(&std::process::id().to_le_bytes());
    h.update(&WRITER_NONCE.fetch_add(1, Ordering::Relaxed).to_le_bytes());
    h.update(&now_wall_ns().to_le_bytes());
    h.finalize().0
}

/// Streaming content-addressed artifact writer (see
/// [`Ledger::artifact_writer`]). Bytes are hashed incrementally and staged
/// as chunk rows under a provisional key inside a writer-owned transaction;
/// `finish` resolves dedupe, rewrites the key to the final content hash, and
/// commits. Dropping without `finish` rolls everything back.
pub struct ArtifactWriter<'a> {
    ledger: &'a Ledger,
    kind: String,
    hasher: Blake3,
    provisional: [u8; 32],
    next_seq: i64,
    buf: Vec<u8>,
    len: u64,
    finished: bool,
}

impl ArtifactWriter<'_> {
    /// Absorb more bytes.
    ///
    /// # Errors
    /// Engine errors while flushing full chunks.
    pub fn write(&mut self, data: &[u8]) -> Result<(), LedgerError> {
        self.hasher.update(data);
        self.len += data.len() as u64;
        self.buf.extend_from_slice(data);
        while self.buf.len() > STORAGE_CHUNK_LEN {
            let rest = self.buf.split_off(STORAGE_CHUNK_LEN);
            let full = std::mem::replace(&mut self.buf, rest);
            self.flush_chunk(&full)?;
        }
        Ok(())
    }

    fn flush_chunk(&mut self, chunk: &[u8]) -> Result<(), LedgerError> {
        self.ledger
            .conn
            .prepare("INSERT INTO artifact_chunks(hash, seq, bytes) VALUES (?1, ?2, ?3)")
            .map_err(|e| sql_err("stream chunk prepare", &e))?
            .execute_with_params(&[
                blob_param(&self.provisional),
                SqliteValue::Integer(self.next_seq),
                blob_param(chunk),
            ])
            .map_err(|e| sql_err("stream chunk insert", &e))?;
        self.next_seq += 1;
        Ok(())
    }

    /// Finalize: dedupe, promote staged chunks to the content hash, commit.
    ///
    /// # Errors
    /// Engine errors; on error the transaction is rolled back and nothing
    /// is stored.
    pub fn finish(mut self, meta: Option<&str>) -> Result<PutReceipt, LedgerError> {
        let result = self.finish_inner(meta);
        self.finished = true;
        if result.is_err() {
            let _ = self.ledger.rollback();
        }
        result
    }

    fn finish_inner(&mut self, meta: Option<&str>) -> Result<PutReceipt, LedgerError> {
        self.ledger.validate_artifact_inputs(&self.kind, meta)?;
        let h = self.hasher.finalize();
        if let Some(info) = self.ledger.artifact_info(&h)? {
            // Identical content already stored: keep theirs — if the
            // bytes still match their identity and the envelope agrees
            // (gp3.19) — then discard staging.
            self.ledger
                .verify_dedupe_candidate(&h, &info, &self.kind, meta)?;
            self.discard_staging()?;
            self.ledger.commit()?;
            return Ok(PutReceipt {
                hash: h,
                len: info.len,
                deduped: true,
                chunked: info.chunk_count > 0,
            });
        }
        if self.next_seq == 0 {
            // Everything fit in the buffer: store inline.
            let buf = std::mem::take(&mut self.buf);
            let receipt = self
                .ledger
                .insert_inline_artifact(&h, &self.kind, &buf, meta)?;
            self.ledger.commit()?;
            return Ok(receipt);
        }
        if !self.buf.is_empty() {
            let tail = std::mem::take(&mut self.buf);
            self.flush_chunk(&tail)?;
        }
        self.ledger
            .conn
            .prepare("UPDATE artifact_chunks SET hash = ?1 WHERE hash = ?2")
            .map_err(|e| sql_err("stream promote prepare", &e))?
            .execute_with_params(&[blob_param(h.as_bytes()), blob_param(&self.provisional)])
            .map_err(|e| sql_err("stream promote", &e))?;
        self.ledger
            .conn
            .prepare(
                "INSERT INTO artifacts(hash, kind, bytes, len, chunk_count, meta, created_at) \
                 VALUES (?1, ?2, NULL, ?3, ?4, ?5, ?6)",
            )
            .map_err(|e| sql_err("stream artifact prepare", &e))?
            .execute_with_params(&[
                blob_param(h.as_bytes()),
                text_param(&self.kind),
                SqliteValue::Integer(int_from_u64(self.len)),
                SqliteValue::Integer(self.next_seq),
                opt_text_param(meta),
                SqliteValue::Integer(now_wall_ns()),
            ])
            .map_err(|e| sql_err("stream artifact insert", &e))?;
        self.ledger.commit()?;
        Ok(PutReceipt {
            hash: h,
            len: self.len,
            deduped: false,
            chunked: true,
        })
    }

    fn discard_staging(&self) -> Result<(), LedgerError> {
        self.ledger
            .conn
            .prepare("DELETE FROM artifact_chunks WHERE hash = ?1")
            .map_err(|e| sql_err("stream discard prepare", &e))?
            .execute_with_params(&[blob_param(&self.provisional)])
            .map_err(|e| sql_err("stream discard", &e))?;
        Ok(())
    }
}

impl Drop for ArtifactWriter<'_> {
    fn drop(&mut self) {
        if !self.finished {
            let _ = self.ledger.rollback();
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn mem() -> Ledger {
        Ledger::open(":memory:").expect("open :memory:")
    }

    fn v7_ledger_with_claim(corrupt_before_migration: bool) -> (Ledger, ContentHash) {
        let conn = Connection::open(":memory:").expect("open v7 fixture");
        conn.query("PRAGMA foreign_keys=ON")
            .expect("enable v7 fixture foreign keys");
        for batch in schema::MIGRATIONS.iter().take(7) {
            for ddl in *batch {
                conn.execute(ddl).expect("apply v7 fixture DDL");
            }
        }
        let instance_bytes = [
            0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x46, 0x77, 0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd,
            0xee, 0xff,
        ];
        conn.prepare("INSERT INTO ledger_identity(singleton, instance_id) VALUES (1, ?1)")
            .expect("prepare v7 fixture identity")
            .execute_with_params(&[blob_param(&instance_bytes)])
            .expect("insert v7 fixture identity");

        let authority = ContentHash([0x41; 32]);
        let governor_hash = ContentHash([0x42; 32]);
        let session_open_hash = ContentHash([0x43; 32]);
        let payload = b"v7-migration-claim".to_vec();
        let source = SessionMutationClaimIdentitySource {
            authority,
            ledger_instance_id: instance_bytes,
            governor_hash,
            session_open_hash,
            registry_schema_version: 1,
            kind: b"meter-report".to_vec(),
            session: 71,
            ledger_scope: b"v7-migration".to_vec(),
            generation: 2,
            causal_ordinal: Some(1),
            payload: payload.clone(),
        };
        let payload_hash = hash_bytes(&payload);
        let claim_hash = ledger_session_mutation_claim_identity(&source);
        conn.prepare(
            "INSERT INTO session_claims( \
                authority, ledger_instance_id, governor_hash, session_open_hash, \
                registry_schema_version, kind, session, ledger_scope, generation, \
                causal_ordinal, payload, payload_hash, claim_hash, created_at \
             ) VALUES (?1, ?2, ?3, ?4, 1, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, 1)",
        )
        .expect("prepare v7 fixture claim")
        .execute_with_params(&[
            blob_param(authority.as_bytes()),
            blob_param(&instance_bytes),
            blob_param(governor_hash.as_bytes()),
            blob_param(session_open_hash.as_bytes()),
            text_param("meter-report"),
            blob_param(&71_u64.to_be_bytes()),
            text_param("v7-migration"),
            blob_param(&2_u64.to_be_bytes()),
            blob_param(&1_u64.to_be_bytes()),
            blob_param(&payload),
            blob_param(payload_hash.as_bytes()),
            blob_param(claim_hash.as_bytes()),
        ])
        .expect("insert v7 fixture claim");
        if corrupt_before_migration {
            conn.execute("DROP TRIGGER trg_session_claims_immutable_update")
                .expect("drop v7 claim update guard for corruption fixture");
            conn.execute("UPDATE session_claims SET kind = 'meter-corrupt'")
                .expect("inject valid-looking v7 semantic corruption");
            let update_guard = schema::V6
                .iter()
                .find(|ddl| {
                    ddl.contains("CREATE TRIGGER IF NOT EXISTS trg_session_claims_immutable_update")
                })
                .expect("shipped v6 claim update guard");
            conn.execute(update_guard)
                .expect("restore exact v7 claim update guard");
        }
        conn.execute("PRAGMA user_version = 7")
            .expect("mark v7 fixture schema");
        (
            Ledger {
                conn,
                path: ":memory:".to_string(),
                instance_id: LedgerInstanceId(instance_bytes),
                read_queries: core::cell::Cell::new(0),
            },
            authority,
        )
    }

    const FX: FiveExplicits<'static> = FiveExplicits {
        seed: &[0x5E, 0xED, 0x00, 0x01],
        versions: r#"{"constellation":"f92683cc4572a198"}"#,
        budget: r#"{"wall_s":10}"#,
        capability: r#"{"ops":["test.*"]}"#,
    };

    fn identity_test_hash(byte: u8) -> ContentHash {
        ContentHash([byte; 32])
    }

    fn assert_identity_moves<T: core::fmt::Debug + PartialEq>(field: &str, base: &T, moved: &T) {
        assert_ne!(base, moved, "semantic identity field {field} did not move");
    }

    #[test]
    fn ledger_semantic_identity_versions_fail_closed() {
        let schemas = [
            (
                PHYSICAL_INSTANCE_IDENTITY_VERSION,
                PHYSICAL_INSTANCE_IDENTITY_DOMAIN,
            ),
            (
                ARTIFACT_CONTENT_IDENTITY_VERSION,
                ARTIFACT_CONTENT_IDENTITY_DOMAIN,
            ),
            (
                SESSION_MUTATION_CLAIM_IDENTITY_VERSION,
                SESSION_MUTATION_CLAIM_IDENTITY_DOMAIN,
            ),
            (
                SESSION_TERMINAL_EVENTS_IDENTITY_VERSION,
                SESSION_TERMINAL_EVENTS_IDENTITY_DOMAIN,
            ),
            (
                SESSION_FLUSH_BATCH_IDENTITY_VERSION,
                SESSION_FLUSH_BATCH_IDENTITY_DOMAIN,
            ),
            (
                SOURCE_ORIGIN_REQUEST_IDENTITY_VERSION,
                SOURCE_ORIGIN_REQUEST_IDENTITY_DOMAIN,
            ),
            (
                DERIVED_COLOR_WAIVER_SUBJECT_IDENTITY_VERSION,
                DERIVED_COLOR_WAIVER_SUBJECT_IDENTITY_DOMAIN,
            ),
            (
                SOURCE_COLOR_WAIVER_SUBJECT_IDENTITY_VERSION,
                SOURCE_COLOR_WAIVER_SUBJECT_IDENTITY_DOMAIN,
            ),
            (COLOR_NODE_IDENTITY_VERSION, COLOR_NODE_IDENTITY_DOMAIN),
            (
                COLOR_ADMISSION_POLICY_IDENTITY_VERSION,
                COLOR_ADMISSION_POLICY_IDENTITY_DOMAIN,
            ),
            (
                VCS_LEDGER_LINEAGE_IDENTITY_VERSION,
                VCS_LEDGER_LINEAGE_IDENTITY_DOMAIN,
            ),
            (
                VCS_COMMIT_LEAF_IDENTITY_VERSION,
                VCS_COMMIT_LEAF_IDENTITY_DOMAIN,
            ),
            (
                VCS_COMMIT_ROOT_IDENTITY_VERSION,
                VCS_COMMIT_ROOT_IDENTITY_DOMAIN,
            ),
            (
                VCS_COMMIT_ENVELOPE_IDENTITY_VERSION,
                VCS_COMMIT_ENVELOPE_IDENTITY_DOMAIN,
            ),
        ];
        for (version, domain) in schemas {
            assert!(identity_schema_is_current(version, domain, version, domain));
            assert!(!identity_schema_is_current(
                version.saturating_sub(1),
                domain,
                version,
                domain
            ));
            assert!(!identity_schema_is_current(
                version + 1,
                domain,
                version,
                domain
            ));
            assert!(!identity_schema_is_current(
                version,
                "org.frankensim.foreign.v1",
                version,
                domain
            ));
        }

        let request = SourceOriginRequestIdentitySource {
            node_name: b"node".to_vec(),
            claimed_color: b"color".to_vec(),
            origin: b"origin".to_vec(),
        };
        let current = ledger_source_origin_request_identity(&request);
        assert!(ledger_source_origin_request_transport_guard(&current));
        let stale = ledger_source_origin_request_identity_with_schema(
            &request,
            0,
            SOURCE_ORIGIN_REQUEST_PREIMAGE_DOMAIN,
        );
        assert!(!ledger_source_origin_request_transport_guard(&stale));

        let derived = derived_waiver_fixture();
        let current = ledger_derived_color_waiver_subject_identity(&derived);
        assert!(ledger_derived_color_waiver_subject_transport_guard(
            &current
        ));
        let stale = ledger_derived_color_waiver_subject_identity_with_schema(
            &derived,
            2,
            COLOR_WAIVER_SUBJECT_PREIMAGE_DOMAIN,
        );
        assert!(!ledger_derived_color_waiver_subject_transport_guard(&stale));

        let source = source_waiver_fixture();
        let current = ledger_source_color_waiver_subject_identity(&source);
        assert!(ledger_source_color_waiver_subject_transport_guard(&current));
        let stale = ledger_source_color_waiver_subject_identity_with_schema(
            &source,
            3,
            COLOR_WAIVER_SUBJECT_PREIMAGE_DOMAIN,
        );
        assert!(!ledger_source_color_waiver_subject_transport_guard(&stale));
    }

    #[test]
    fn physical_instance_identity_fields_move_independently() {
        let source = PhysicalInstanceIdentitySource { uuid: [1; 16] };
        let base = ledger_physical_instance_identity(&source);
        let mut changed = source.clone();
        changed.uuid[15] ^= 1;
        assert_identity_moves(
            "uuid-bytes",
            &base,
            &ledger_physical_instance_identity(&changed),
        );
    }

    #[test]
    fn physical_instance_excluded_fields_do_not_move_identity() {
        let source = PhysicalInstanceIdentitySource { uuid: [2; 16] };
        let base = ledger_physical_instance_identity(&source);
        let envelope_a = ("/old/path", "/alias/a", 0x1000_usize, 1_u64);
        let envelope_b = ("/new/path", "/alias/b", 0x2000_usize, 9_u64);
        assert_ne!(envelope_a, envelope_b);
        assert_eq!(base, ledger_physical_instance_identity(&source));
    }

    #[test]
    fn artifact_content_identity_fields_move_independently() {
        let source = ArtifactContentIdentitySource {
            content: b"artifact".to_vec(),
        };
        let base = ledger_artifact_content_identity(&source);
        let mut changed = source.clone();
        changed.content.push(0);
        assert_identity_moves(
            "content-bytes",
            &base,
            &ledger_artifact_content_identity(&changed),
        );
    }

    #[test]
    fn artifact_content_excluded_fields_do_not_move_identity() {
        let source = ArtifactContentIdentitySource {
            content: b"same-content".to_vec(),
        };
        let base = ledger_artifact_content_identity(&source);
        let envelope_a = ("mesh", "{}", 1_i64, vec![4_usize, 8]);
        let envelope_b = ("field", "{\"unit\":\"m\"}", 2_i64, vec![3_usize, 9]);
        assert_ne!(envelope_a, envelope_b);
        assert_eq!(base, ledger_artifact_content_identity(&source));
    }

    fn session_claim_fixture() -> SessionMutationClaimIdentitySource {
        SessionMutationClaimIdentitySource {
            authority: identity_test_hash(1),
            ledger_instance_id: [2; 16],
            governor_hash: identity_test_hash(3),
            session_open_hash: identity_test_hash(4),
            registry_schema_version: 1,
            kind: b"submission".to_vec(),
            session: 5,
            ledger_scope: b"scope".to_vec(),
            generation: 6,
            causal_ordinal: Some(7),
            payload: b"payload".to_vec(),
        }
    }

    #[test]
    fn session_mutation_claim_identity_fields_move_independently() {
        let source = session_claim_fixture();
        let base = ledger_session_mutation_claim_identity(&source);
        let foreign_domain = ledger_session_mutation_claim_identity_with_domain(
            &source,
            b"org.frankensim.fs-ledger.session-mutation-claim.w1\0",
        );
        assert_identity_moves("identity-domain", &base, &foreign_domain);

        let mut changed = source.clone();
        changed.authority = identity_test_hash(10);
        assert_identity_moves(
            "authority",
            &base,
            &ledger_session_mutation_claim_identity(&changed),
        );
        changed = source.clone();
        changed.ledger_instance_id[0] ^= 1;
        assert_identity_moves(
            "ledger-instance-id",
            &base,
            &ledger_session_mutation_claim_identity(&changed),
        );
        changed = source.clone();
        changed.governor_hash = identity_test_hash(11);
        assert_identity_moves(
            "governor-hash",
            &base,
            &ledger_session_mutation_claim_identity(&changed),
        );
        changed = source.clone();
        changed.session_open_hash = identity_test_hash(12);
        assert_identity_moves(
            "session-open-hash",
            &base,
            &ledger_session_mutation_claim_identity(&changed),
        );
        changed = source.clone();
        changed.registry_schema_version += 1;
        assert_identity_moves(
            "registry-schema-version",
            &base,
            &ledger_session_mutation_claim_identity(&changed),
        );
        changed = source.clone();
        changed.kind.push(b'x');
        for field in ["kind-byte-count", "kind-bytes"] {
            assert_identity_moves(
                field,
                &base,
                &ledger_session_mutation_claim_identity(&changed),
            );
        }
        changed = source.clone();
        changed.session += 1;
        assert_identity_moves(
            "session",
            &base,
            &ledger_session_mutation_claim_identity(&changed),
        );
        changed = source.clone();
        changed.ledger_scope.push(b'x');
        for field in ["ledger-scope-byte-count", "ledger-scope-bytes"] {
            assert_identity_moves(
                field,
                &base,
                &ledger_session_mutation_claim_identity(&changed),
            );
        }
        changed = source.clone();
        changed.generation += 1;
        assert_identity_moves(
            "generation",
            &base,
            &ledger_session_mutation_claim_identity(&changed),
        );
        changed = source.clone();
        changed.causal_ordinal = None;
        assert_identity_moves(
            "causal-ordinal-presence",
            &base,
            &ledger_session_mutation_claim_identity(&changed),
        );
        changed = source.clone();
        changed.causal_ordinal = Some(8);
        assert_identity_moves(
            "causal-ordinal-value",
            &base,
            &ledger_session_mutation_claim_identity(&changed),
        );
        changed = source;
        changed.payload.push(0);
        assert_identity_moves(
            "payload-bytes-via-blake3",
            &base,
            &ledger_session_mutation_claim_identity(&changed),
        );
    }

    #[test]
    fn session_mutation_claim_excluded_fields_do_not_move_identity() {
        let source = session_claim_fixture();
        let base = ledger_session_mutation_claim_identity(&source);
        let envelope_a = (1_i64, 10_i64, identity_test_hash(1));
        let envelope_b = (2_i64, 20_i64, identity_test_hash(2));
        assert_ne!(envelope_a, envelope_b);
        assert_eq!(base, ledger_session_mutation_claim_identity(&source));
    }

    fn session_events_fixture() -> Vec<SessionTerminalEventIdentitySource> {
        vec![
            SessionTerminalEventIdentitySource {
                session: 1_u64.to_be_bytes().to_vec(),
                timestamp: 10,
                kind: b"first".to_vec(),
                payload: Some(b"{\"v\":1}".to_vec()),
            },
            SessionTerminalEventIdentitySource {
                session: 1_u64.to_be_bytes().to_vec(),
                timestamp: 11,
                kind: b"second".to_vec(),
                payload: None,
            },
        ]
    }

    #[test]
    fn session_terminal_events_identity_fields_move_independently() {
        let events = session_events_fixture();
        let base = ledger_session_terminal_events_identity(&events);
        assert_identity_moves(
            "identity-domain",
            &base,
            &ledger_session_terminal_events_identity_with_schema(
                &events,
                events.len(),
                b"org.frankensim.fs-ledger.session-terminal-events.w2\0",
            ),
        );
        assert_identity_moves(
            "event-count",
            &base,
            &ledger_session_terminal_events_identity_with_schema(
                &events,
                events.len() + 1,
                b"org.frankensim.fs-ledger.session-terminal-events.v2\0",
            ),
        );
        let mut changed = events.clone();
        changed.swap(0, 1);
        assert_identity_moves(
            "event-order",
            &base,
            &ledger_session_terminal_events_identity(&changed),
        );
        changed = events.clone();
        changed[0].session.push(0);
        for field in ["session-byte-count", "session-bytes"] {
            assert_identity_moves(
                field,
                &base,
                &ledger_session_terminal_events_identity(&changed),
            );
        }
        changed = events.clone();
        changed[0].timestamp += 1;
        assert_identity_moves(
            "timestamp",
            &base,
            &ledger_session_terminal_events_identity(&changed),
        );
        changed = events.clone();
        changed[0].kind.push(b'x');
        for field in ["kind-byte-count", "kind-bytes"] {
            assert_identity_moves(
                field,
                &base,
                &ledger_session_terminal_events_identity(&changed),
            );
        }
        changed = events.clone();
        changed[0].payload = None;
        assert_identity_moves(
            "payload-presence",
            &base,
            &ledger_session_terminal_events_identity(&changed),
        );
        changed = events;
        changed[0]
            .payload
            .as_mut()
            .expect("fixture payload")
            .push(b' ');
        for field in ["payload-byte-count", "payload-bytes"] {
            assert_identity_moves(
                field,
                &base,
                &ledger_session_terminal_events_identity(&changed),
            );
        }
    }

    #[test]
    fn session_terminal_events_excluded_fields_do_not_move_identity() {
        let events = session_events_fixture();
        let base = ledger_session_terminal_events_identity(&events);
        let storage_a = (1_i64, 2_i64, 3_i64);
        let storage_b = (4_i64, 5_i64, 6_i64);
        assert_ne!(storage_a, storage_b);
        assert_eq!(base, ledger_session_terminal_events_identity(&events));
    }

    fn session_batch_fixture() -> SessionFlushBatchIdentitySource {
        SessionFlushBatchIdentitySource {
            ledger_instance_id: [1; 16],
            registry_schema_version: 1,
            terminals: vec![
                SessionFlushTerminalIdentitySource {
                    authority: identity_test_hash(2),
                    claim_hash: identity_test_hash(3),
                    receipt_hash: identity_test_hash(4),
                    event_count: 1,
                    events_hash: identity_test_hash(5),
                    encoded_bytes: 100,
                },
                SessionFlushTerminalIdentitySource {
                    authority: identity_test_hash(6),
                    claim_hash: identity_test_hash(7),
                    receipt_hash: identity_test_hash(8),
                    event_count: 2,
                    events_hash: identity_test_hash(9),
                    encoded_bytes: 200,
                },
            ],
        }
    }

    #[test]
    fn session_flush_batch_identity_fields_move_independently() {
        let source = session_batch_fixture();
        let base = ledger_session_flush_batch_identity(&source);
        assert_identity_moves(
            "identity-domain",
            &base,
            &ledger_session_flush_batch_identity_with_schema(
                &source,
                source.terminals.len(),
                b"org.frankensim.fs-ledger.session-flush-batch.w2\0",
            ),
        );
        let mut changed = source.clone();
        changed.ledger_instance_id[0] ^= 1;
        assert_identity_moves(
            "ledger-instance-id",
            &base,
            &ledger_session_flush_batch_identity(&changed),
        );
        changed = source.clone();
        changed.registry_schema_version += 1;
        assert_identity_moves(
            "registry-schema-version",
            &base,
            &ledger_session_flush_batch_identity(&changed),
        );
        assert_identity_moves(
            "terminal-count",
            &base,
            &ledger_session_flush_batch_identity_with_schema(
                &source,
                source.terminals.len() + 1,
                b"org.frankensim.fs-ledger.session-flush-batch.v2\0",
            ),
        );
        changed = source.clone();
        changed.terminals.swap(0, 1);
        assert_identity_moves(
            "terminal-order",
            &base,
            &ledger_session_flush_batch_identity(&changed),
        );
        let field_mutations: [(&str, fn(&mut SessionFlushTerminalIdentitySource)); 6] = [
            ("authority", |row| row.authority = identity_test_hash(20)),
            ("claim-hash", |row| row.claim_hash = identity_test_hash(21)),
            ("receipt-hash", |row| {
                row.receipt_hash = identity_test_hash(22);
            }),
            ("event-count", |row| row.event_count += 1),
            ("events-hash", |row| {
                row.events_hash = identity_test_hash(23)
            }),
            ("encoded-byte-count", |row| row.encoded_bytes += 1),
        ];
        for (field, mutate) in field_mutations {
            changed = source.clone();
            mutate(&mut changed.terminals[0]);
            assert_identity_moves(field, &base, &ledger_session_flush_batch_identity(&changed));
        }
    }

    #[test]
    fn session_flush_batch_excluded_fields_do_not_move_identity() {
        let source = session_batch_fixture();
        let base = ledger_session_flush_batch_identity(&source);
        let storage_a = (1_i64, 2_i64, 3_i64, identity_test_hash(4));
        let storage_b = (5_i64, 6_i64, 7_i64, identity_test_hash(8));
        assert_ne!(storage_a, storage_b);
        assert_eq!(base, ledger_session_flush_batch_identity(&source));
    }

    #[test]
    fn source_origin_request_identity_fields_move_independently() {
        let source = SourceOriginRequestIdentitySource {
            node_name: b"node".to_vec(),
            claimed_color: b"exact-color".to_vec(),
            origin: b"typed-origin".to_vec(),
        };
        let base = ledger_source_origin_request_identity(&source);
        assert_identity_moves(
            "transport-version",
            &base,
            &ledger_source_origin_request_identity_with_schema(
                &source,
                2,
                SOURCE_ORIGIN_REQUEST_PREIMAGE_DOMAIN,
            ),
        );
        assert_identity_moves(
            "preimage-domain",
            &base,
            &ledger_source_origin_request_identity_with_schema(
                &source,
                1,
                b"frankensim/fs-ledger/source-origin-requesu",
            ),
        );
        let mut corrupt_count = base.clone();
        corrupt_count[1] ^= 1;
        assert_identity_moves("domain-byte-count", &base, &corrupt_count);

        let mut changed = source.clone();
        changed.node_name.push(b'x');
        for field in ["node-name-byte-count", "node-name"] {
            assert_identity_moves(
                field,
                &base,
                &ledger_source_origin_request_identity(&changed),
            );
        }
        changed = source.clone();
        changed.claimed_color.push(0);
        for field in ["claimed-color-byte-count", "claimed-color-canonical-bytes"] {
            assert_identity_moves(
                field,
                &base,
                &ledger_source_origin_request_identity(&changed),
            );
        }
        changed = source;
        changed.origin.push(0);
        assert_identity_moves(
            "typed-origin-canonical-bytes",
            &base,
            &ledger_source_origin_request_identity(&changed),
        );
    }

    #[test]
    fn source_origin_request_excluded_fields_do_not_move_identity() {
        let source = SourceOriginRequestIdentitySource {
            node_name: b"node".to_vec(),
            claimed_color: b"color".to_vec(),
            origin: b"origin".to_vec(),
        };
        let base = ledger_source_origin_request_identity(&source);
        let callback_a = (identity_test_hash(1), 1_usize);
        let callback_b = (identity_test_hash(2), 9_usize);
        assert_ne!(callback_a, callback_b);
        assert_eq!(base, ledger_source_origin_request_identity(&source));
    }

    fn derived_waiver_fixture() -> DerivedColorWaiverSubjectIdentitySource {
        DerivedColorWaiverSubjectIdentitySource {
            operation_tag: 1,
            key_id: b"key".to_vec(),
            scope: b"color-upgrade".to_vec(),
            node_name: b"node".to_vec(),
            claimed_color: b"color".to_vec(),
            annotation_id: b"waiver".to_vec(),
            annotation_signer: b"signer".to_vec(),
            annotation_reason: b"reason".to_vec(),
            parent_hashes: vec![identity_test_hash(1), identity_test_hash(2)],
            expires_day: 100,
            signature: b"signature".to_vec(),
        }
    }

    fn source_waiver_fixture() -> SourceColorWaiverSubjectIdentitySource {
        SourceColorWaiverSubjectIdentitySource {
            key_id: b"key".to_vec(),
            scope: b"source-color".to_vec(),
            node_name: b"node".to_vec(),
            claimed_color: b"color".to_vec(),
            annotation_id: b"waiver".to_vec(),
            annotation_signer: b"signer".to_vec(),
            annotation_reason: b"reason".to_vec(),
            parent_hashes: vec![identity_test_hash(1), identity_test_hash(2)],
            expires_day: 100,
            signature: b"signature".to_vec(),
        }
    }

    #[test]
    fn derived_color_waiver_subject_identity_fields_move_independently() {
        let source = derived_waiver_fixture();
        let base = ledger_derived_color_waiver_subject_identity(&source);
        assert_identity_moves(
            "transport-version",
            &base,
            &ledger_derived_color_waiver_subject_identity_with_schema(
                &source,
                2,
                COLOR_WAIVER_SUBJECT_PREIMAGE_DOMAIN,
            ),
        );
        assert_identity_moves(
            "preimage-domain",
            &base,
            &ledger_derived_color_waiver_subject_identity_with_schema(
                &source,
                3,
                b"frankensim/fs-ledger/color-waiveq",
            ),
        );
        let mut corrupt_count = base.clone();
        corrupt_count[1] ^= 1;
        assert_identity_moves("domain-byte-count", &base, &corrupt_count);

        let mut changed = source.clone();
        changed.operation_tag ^= 1;
        assert_identity_moves(
            "operation-tag",
            &base,
            &ledger_derived_color_waiver_subject_identity(&changed),
        );
        let text_mutations: [(&[&str], fn(&mut DerivedColorWaiverSubjectIdentitySource)); 7] = [
            (&["key-id-byte-count", "key-id"], |value| {
                value.key_id.push(b'x')
            }),
            (&["scope-byte-count", "scope"], |value| {
                value.scope.push(b'x')
            }),
            (&["node-name-byte-count", "node-name"], |value| {
                value.node_name.push(b'x')
            }),
            (
                &["claimed-color-byte-count", "claimed-color-canonical-bytes"],
                |value| value.claimed_color.push(0),
            ),
            (&["annotation-id-byte-count", "annotation-id"], |value| {
                value.annotation_id.push(b'x')
            }),
            (
                &["annotation-signer-byte-count", "annotation-signer"],
                |value| value.annotation_signer.push(b'x'),
            ),
            (
                &["annotation-reason-byte-count", "annotation-reason"],
                |value| value.annotation_reason.push(b'x'),
            ),
        ];
        for (fields, mutate) in text_mutations {
            changed = source.clone();
            mutate(&mut changed);
            let moved = ledger_derived_color_waiver_subject_identity(&changed);
            for field in fields {
                assert_identity_moves(field, &base, &moved);
            }
        }
        changed = source.clone();
        changed.parent_hashes.push(identity_test_hash(3));
        assert_identity_moves(
            "parent-count",
            &base,
            &ledger_derived_color_waiver_subject_identity(&changed),
        );
        changed = source.clone();
        changed.parent_hashes.swap(0, 1);
        assert_identity_moves(
            "parent-order",
            &base,
            &ledger_derived_color_waiver_subject_identity(&changed),
        );
        changed = source.clone();
        changed.parent_hashes[0] = identity_test_hash(4);
        assert_identity_moves(
            "parent-hashes",
            &base,
            &ledger_derived_color_waiver_subject_identity(&changed),
        );
        changed = source;
        changed.expires_day += 1;
        assert_identity_moves(
            "expires-day",
            &base,
            &ledger_derived_color_waiver_subject_identity(&changed),
        );
    }

    #[test]
    fn derived_color_waiver_subject_excluded_fields_do_not_move_identity() {
        let source = derived_waiver_fixture();
        let base = ledger_derived_color_waiver_subject_identity(&source);
        let mut changed = source.clone();
        changed.signature.push(0);
        assert_eq!(
            base,
            ledger_derived_color_waiver_subject_identity(&changed),
            "a signature is not part of its own signing subject"
        );
        let admission_a = (1_u32, identity_test_hash(1));
        let admission_b = (2_u32, identity_test_hash(2));
        assert_ne!(admission_a, admission_b);
        assert_eq!(base, ledger_derived_color_waiver_subject_identity(&source));
    }

    #[test]
    fn source_color_waiver_subject_identity_fields_move_independently() {
        let source = source_waiver_fixture();
        let base = ledger_source_color_waiver_subject_identity(&source);
        assert_identity_moves(
            "transport-version",
            &base,
            &ledger_source_color_waiver_subject_identity_with_schema(
                &source,
                3,
                COLOR_WAIVER_SUBJECT_PREIMAGE_DOMAIN,
            ),
        );
        assert_identity_moves(
            "preimage-domain",
            &base,
            &ledger_source_color_waiver_subject_identity_with_schema(
                &source,
                4,
                b"frankensim/fs-ledger/color-waiveq",
            ),
        );
        let mut corrupt_count = base.clone();
        corrupt_count[1] ^= 1;
        assert_identity_moves("domain-byte-count", &base, &corrupt_count);
        let sentinel = 1 + 8 + COLOR_WAIVER_SUBJECT_PREIMAGE_DOMAIN.len();
        let mut changed_sentinel = base.clone();
        changed_sentinel[sentinel] ^= 1;
        assert_identity_moves("source-operation-sentinel", &base, &changed_sentinel);

        let text_mutations: [(&[&str], fn(&mut SourceColorWaiverSubjectIdentitySource)); 7] = [
            (&["key-id-byte-count", "key-id"], |value| {
                value.key_id.push(b'x')
            }),
            (&["scope-byte-count", "scope"], |value| {
                value.scope.push(b'x')
            }),
            (&["node-name-byte-count", "node-name"], |value| {
                value.node_name.push(b'x')
            }),
            (
                &["claimed-color-byte-count", "claimed-color-canonical-bytes"],
                |value| value.claimed_color.push(0),
            ),
            (&["annotation-id-byte-count", "annotation-id"], |value| {
                value.annotation_id.push(b'x')
            }),
            (
                &["annotation-signer-byte-count", "annotation-signer"],
                |value| value.annotation_signer.push(b'x'),
            ),
            (
                &["annotation-reason-byte-count", "annotation-reason"],
                |value| value.annotation_reason.push(b'x'),
            ),
        ];
        for (fields, mutate) in text_mutations {
            let mut changed = source.clone();
            mutate(&mut changed);
            let moved = ledger_source_color_waiver_subject_identity(&changed);
            for field in fields {
                assert_identity_moves(field, &base, &moved);
            }
        }
        let mut changed = source.clone();
        changed.parent_hashes.push(identity_test_hash(3));
        assert_identity_moves(
            "parent-count",
            &base,
            &ledger_source_color_waiver_subject_identity(&changed),
        );
        changed = source.clone();
        changed.parent_hashes.swap(0, 1);
        assert_identity_moves(
            "parent-order",
            &base,
            &ledger_source_color_waiver_subject_identity(&changed),
        );
        changed = source.clone();
        changed.parent_hashes[0] = identity_test_hash(4);
        assert_identity_moves(
            "parent-hashes",
            &base,
            &ledger_source_color_waiver_subject_identity(&changed),
        );
        changed = source;
        changed.expires_day += 1;
        assert_identity_moves(
            "expires-day",
            &base,
            &ledger_source_color_waiver_subject_identity(&changed),
        );
    }

    #[test]
    fn source_color_waiver_subject_excluded_fields_do_not_move_identity() {
        let source = source_waiver_fixture();
        let base = ledger_source_color_waiver_subject_identity(&source);
        let mut changed = source.clone();
        changed.signature.push(0);
        assert_eq!(
            base,
            ledger_source_color_waiver_subject_identity(&changed),
            "a signature is not part of its own signing subject"
        );
        let admission_a = (1_u32, identity_test_hash(1));
        let admission_b = (2_u32, identity_test_hash(2));
        assert_ne!(admission_a, admission_b);
        assert_eq!(base, ledger_source_color_waiver_subject_identity(&source));
    }

    fn color_node_fixture() -> ColorNodeIdentitySource {
        ColorNodeIdentitySource {
            node_id: 10,
            operation_tag: Some(1),
            name: b"node".to_vec(),
            color: b"canonical-color".to_vec(),
            parent_local_ids: vec![1, 2],
            parent_hashes: vec![identity_test_hash(1), identity_test_hash(2)],
            demotions: vec![b"demotion-a".to_vec(), b"demotion-b".to_vec()],
            origin: Some(b"origin".to_vec()),
            origin_policy_fingerprint: Some(identity_test_hash(3)),
            waiver_dependencies: vec![b"dependency-a".to_vec(), b"dependency-b".to_vec()],
            waiver: Some(b"waiver".to_vec()),
            grant_payload: Some(b"grant-payload".to_vec()),
            grant_signature: Some(b"grant-signature".to_vec()),
            waiver_policy_fingerprint: Some(identity_test_hash(4)),
            waiver_admission_day: Some(20),
            stored_hash: identity_test_hash(5),
        }
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn color_node_identity_fields_move_independently() {
        let source = color_node_fixture();
        let base = ledger_color_node_identity(&source);
        assert_identity_moves(
            "transport-version",
            &base,
            &ledger_color_node_identity_with_schema(
                &source,
                COLOR_NODE_IDENTITY_VERSION as u8 + 1,
                COLOR_NODE_PREIMAGE_DOMAIN,
            ),
        );
        let domain_moved = ledger_color_node_identity_with_schema(
            &source,
            COLOR_NODE_IDENTITY_VERSION as u8,
            b"frankensim/fs-ledger/color-node/w2",
        );
        for field in ["domain-byte-count", "preimage-domain"] {
            assert_identity_moves(field, &base, &domain_moved);
        }
        let mut changed = source.clone();
        changed.operation_tag = None;
        assert_identity_moves(
            "operation-presence",
            &base,
            &ledger_color_node_identity(&changed),
        );
        changed = source.clone();
        changed.operation_tag = Some(2);
        assert_identity_moves(
            "operation-tag",
            &base,
            &ledger_color_node_identity(&changed),
        );
        changed = source.clone();
        changed.name.push(b'x');
        for field in ["node-name-byte-count", "node-name"] {
            assert_identity_moves(field, &base, &ledger_color_node_identity(&changed));
        }
        changed = source.clone();
        changed.color.push(0);
        for field in ["color-byte-count", "color-canonical-bytes"] {
            assert_identity_moves(field, &base, &ledger_color_node_identity(&changed));
        }
        changed = source.clone();
        changed.parent_hashes.push(identity_test_hash(6));
        assert_identity_moves("parent-count", &base, &ledger_color_node_identity(&changed));
        changed = source.clone();
        changed.parent_hashes.swap(0, 1);
        assert_identity_moves("parent-order", &base, &ledger_color_node_identity(&changed));
        changed = source.clone();
        changed.parent_hashes[0] = identity_test_hash(7);
        assert_identity_moves(
            "parent-hashes",
            &base,
            &ledger_color_node_identity(&changed),
        );
        changed = source.clone();
        changed.demotions.push(b"demotion-c".to_vec());
        assert_identity_moves(
            "demotion-count",
            &base,
            &ledger_color_node_identity(&changed),
        );
        changed = source.clone();
        changed.demotions.swap(0, 1);
        assert_identity_moves(
            "demotion-order",
            &base,
            &ledger_color_node_identity(&changed),
        );
        changed = source.clone();
        changed.demotions[0].push(0);
        assert_identity_moves(
            "demotion-canonical-bytes",
            &base,
            &ledger_color_node_identity(&changed),
        );
        changed = source.clone();
        changed.origin = None;
        assert_identity_moves(
            "origin-presence",
            &base,
            &ledger_color_node_identity(&changed),
        );
        changed = source.clone();
        changed.origin.as_mut().expect("fixture origin").push(0);
        assert_identity_moves(
            "origin-canonical-bytes",
            &base,
            &ledger_color_node_identity(&changed),
        );
        changed = source.clone();
        changed.origin_policy_fingerprint = None;
        assert_identity_moves(
            "origin-policy-presence",
            &base,
            &ledger_color_node_identity(&changed),
        );
        changed = source.clone();
        changed.origin_policy_fingerprint = Some(identity_test_hash(8));
        assert_identity_moves(
            "origin-policy-fingerprint",
            &base,
            &ledger_color_node_identity(&changed),
        );
        changed = source.clone();
        changed.waiver_dependencies.push(b"dependency-c".to_vec());
        assert_identity_moves(
            "waiver-dependency-count",
            &base,
            &ledger_color_node_identity(&changed),
        );
        changed = source.clone();
        changed.waiver_dependencies.swap(0, 1);
        assert_identity_moves(
            "waiver-dependency-order",
            &base,
            &ledger_color_node_identity(&changed),
        );
        changed = source.clone();
        changed.waiver_dependencies[0].push(0);
        assert_identity_moves(
            "waiver-dependency-canonical-bytes",
            &base,
            &ledger_color_node_identity(&changed),
        );
        changed = source.clone();
        changed.waiver = None;
        assert_identity_moves(
            "waiver-presence",
            &base,
            &ledger_color_node_identity(&changed),
        );
        changed = source.clone();
        changed.waiver.as_mut().expect("fixture waiver").push(0);
        assert_identity_moves(
            "waiver-canonical-bytes",
            &base,
            &ledger_color_node_identity(&changed),
        );
        changed = source.clone();
        changed.grant_payload = None;
        assert_identity_moves(
            "grant-presence",
            &base,
            &ledger_color_node_identity(&changed),
        );
        changed = source.clone();
        changed
            .grant_payload
            .as_mut()
            .expect("fixture grant payload")
            .push(0);
        assert_identity_moves(
            "grant-payload",
            &base,
            &ledger_color_node_identity(&changed),
        );
        changed = source.clone();
        changed
            .grant_signature
            .as_mut()
            .expect("fixture grant signature")
            .push(0);
        assert_identity_moves(
            "grant-signature",
            &base,
            &ledger_color_node_identity(&changed),
        );
        changed = source.clone();
        changed.waiver_policy_fingerprint = None;
        assert_identity_moves(
            "waiver-policy-presence",
            &base,
            &ledger_color_node_identity(&changed),
        );
        changed = source.clone();
        changed.waiver_policy_fingerprint = Some(identity_test_hash(9));
        assert_identity_moves(
            "waiver-policy-fingerprint",
            &base,
            &ledger_color_node_identity(&changed),
        );
        changed = source.clone();
        changed.waiver_admission_day = None;
        assert_identity_moves(
            "waiver-admission-day-presence",
            &base,
            &ledger_color_node_identity(&changed),
        );
        changed = source;
        changed.waiver_admission_day = Some(21);
        assert_identity_moves(
            "waiver-admission-day",
            &base,
            &ledger_color_node_identity(&changed),
        );
    }

    #[test]
    fn color_node_excluded_fields_do_not_move_identity() {
        let source = color_node_fixture();
        let base = ledger_color_node_identity(&source);
        let mut changed = source.clone();
        changed.node_id += 1;
        changed.parent_local_ids = vec![100, 200];
        changed.stored_hash = identity_test_hash(99);
        assert_eq!(base, ledger_color_node_identity(&changed));
        let envelope_a = ("rounded-json-a", 10_i64, 6_u32);
        let envelope_b = ("rounded-json-b", 20_i64, 7_u32);
        assert_ne!(envelope_a, envelope_b);
        assert_eq!(base, ledger_color_node_identity(&source));
    }

    #[test]
    fn color_admission_policy_identity_fields_move_independently() {
        let source = ColorAdmissionPolicyIdentitySource {
            color_write_row_schema_version: 7,
            color_algebra_version: 2,
        };
        let base = ledger_color_admission_policy_identity(&source);
        assert_identity_moves(
            "preimage-domain",
            &base,
            &ledger_color_admission_policy_identity_with_schema(
                &source,
                "fs-ledger/color-admission-policy/w1",
            ),
        );
        let mut changed = source.clone();
        changed.color_write_row_schema_version += 1;
        assert_identity_moves(
            "color-write-row-schema-version",
            &base,
            &ledger_color_admission_policy_identity(&changed),
        );
        changed = source;
        changed.color_algebra_version += 1;
        assert_identity_moves(
            "color-algebra-version",
            &base,
            &ledger_color_admission_policy_identity(&changed),
        );
    }

    #[test]
    fn color_admission_policy_excluded_fields_do_not_move_identity() {
        let source = ColorAdmissionPolicyIdentitySource {
            color_write_row_schema_version: 7,
            color_algebra_version: 2,
        };
        let base = ledger_color_admission_policy_identity(&source);
        let envelope_a = ("build-a", 10_i64);
        let envelope_b = ("build-b", 20_i64);
        assert_ne!(envelope_a, envelope_b);
        assert_eq!(base, ledger_color_admission_policy_identity(&source));
    }

    #[test]
    fn vcs_ledger_lineage_identity_fields_move_independently() {
        let source = VcsLedgerLineageIdentitySource {
            mint_path: b"/ledger/path".to_vec(),
            minted_ns: 10,
        };
        let base = ledger_vcs_ledger_lineage_identity(&source);
        let moved_domain = ledger_vcs_ledger_lineage_identity_with_domain(
            &source,
            b"frankensim.fs-ledger.vcs.ledger-identity.w1",
        );
        for field in ["domain-label-frame", "domain-byte-count", "preimage-domain"] {
            assert_identity_moves(field, &base, &moved_domain);
        }
        let mut changed = source.clone();
        changed.mint_path.push(b'x');
        for field in ["mint-path-byte-count", "mint-path"] {
            assert_identity_moves(field, &base, &ledger_vcs_ledger_lineage_identity(&changed));
        }
        changed = source;
        changed.minted_ns += 1;
        assert_identity_moves(
            "minted-nanoseconds",
            &base,
            &ledger_vcs_ledger_lineage_identity(&changed),
        );
    }

    #[test]
    fn vcs_ledger_lineage_excluded_fields_do_not_move_identity() {
        let source = VcsLedgerLineageIdentitySource {
            mint_path: b"/ledger/path".to_vec(),
            minted_ns: 10,
        };
        let base = ledger_vcs_ledger_lineage_identity(&source);
        let persisted_a = (1_i64, "/copy/a", 20_i64);
        let persisted_b = (9_i64, "/copy/b", 30_i64);
        assert_ne!(persisted_a, persisted_b);
        assert_eq!(base, ledger_vcs_ledger_lineage_identity(&source));
    }

    fn vcs_leaf_fixture() -> VcsCommitLeafIdentitySource {
        VcsCommitLeafIdentitySource {
            ir: b"{\"op\":\"fixture\"}".to_vec(),
            seed: b"seed".to_vec(),
            versions: b"{\"v\":1}".to_vec(),
            budget: b"{\"wall\":1}".to_vec(),
            capability: b"{\"ops\":[\"fixture\"]}".to_vec(),
            outcome: Some(b"ok".to_vec()),
            diagnostic: Some(b"{\"detail\":\"done\"}".to_vec()),
            execution_mode: b"deterministic".to_vec(),
            edges: vec![
                VcsCommitEdgeIdentitySource {
                    role: b"in".to_vec(),
                    artifact_hash: identity_test_hash(1),
                },
                VcsCommitEdgeIdentitySource {
                    role: b"out".to_vec(),
                    artifact_hash: identity_test_hash(2),
                },
            ],
        }
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn vcs_commit_leaf_identity_fields_move_independently() {
        let source = vcs_leaf_fixture();
        let base = ledger_vcs_commit_leaf_identity(&source);
        let moved_domain = ledger_vcs_commit_leaf_identity_with_domain(
            &source,
            b"frankensim.fs-ledger.vcs.commit-leaf.w2",
        );
        for field in ["domain-label-frame", "domain-byte-count", "preimage-domain"] {
            assert_identity_moves(field, &base, &moved_domain);
        }
        let text_mutations: [(&[&str], fn(&mut VcsCommitLeafIdentitySource)); 6] = [
            (&["ir-byte-count", "ir-bytes"], |value| value.ir.push(b' ')),
            (&["seed-byte-count", "seed-bytes"], |value| {
                value.seed.push(0)
            }),
            (&["versions-byte-count", "versions-bytes"], |value| {
                value.versions.push(b' ')
            }),
            (&["budget-byte-count", "budget-bytes"], |value| {
                value.budget.push(b' ')
            }),
            (&["capability-byte-count", "capability-bytes"], |value| {
                value.capability.push(b' ')
            }),
            (&["execution-mode-byte-count", "execution-mode"], |value| {
                value.execution_mode.push(b'x')
            }),
        ];
        for (fields, mutate) in text_mutations {
            let mut changed = source.clone();
            mutate(&mut changed);
            let moved = ledger_vcs_commit_leaf_identity(&changed);
            for field in fields {
                assert_identity_moves(field, &base, &moved);
            }
        }
        let mut changed = source.clone();
        changed.outcome = None;
        assert_identity_moves(
            "outcome-presence",
            &base,
            &ledger_vcs_commit_leaf_identity(&changed),
        );
        changed = source.clone();
        changed
            .outcome
            .as_mut()
            .expect("fixture outcome")
            .push(b'x');
        for field in ["outcome-byte-count", "outcome-bytes"] {
            assert_identity_moves(field, &base, &ledger_vcs_commit_leaf_identity(&changed));
        }
        changed = source.clone();
        changed.diagnostic = None;
        assert_identity_moves(
            "diagnostic-presence",
            &base,
            &ledger_vcs_commit_leaf_identity(&changed),
        );
        changed = source.clone();
        changed
            .diagnostic
            .as_mut()
            .expect("fixture diagnostic")
            .push(b' ');
        for field in ["diagnostic-byte-count", "diagnostic-bytes"] {
            assert_identity_moves(field, &base, &ledger_vcs_commit_leaf_identity(&changed));
        }
        changed = source.clone();
        changed.edges.push(VcsCommitEdgeIdentitySource {
            role: b"out".to_vec(),
            artifact_hash: identity_test_hash(3),
        });
        assert_identity_moves(
            "edge-count",
            &base,
            &ledger_vcs_commit_leaf_identity(&changed),
        );
        changed = source.clone();
        changed.edges.swap(0, 1);
        assert_identity_moves(
            "edge-order",
            &base,
            &ledger_vcs_commit_leaf_identity(&changed),
        );
        changed = source.clone();
        changed.edges[0].role.push(b'x');
        for field in ["edge-role-byte-count", "edge-role"] {
            assert_identity_moves(field, &base, &ledger_vcs_commit_leaf_identity(&changed));
        }
        changed = source;
        changed.edges[0].artifact_hash = identity_test_hash(4);
        assert_identity_moves(
            "artifact-hash",
            &base,
            &ledger_vcs_commit_leaf_identity(&changed),
        );
    }

    #[test]
    fn vcs_commit_leaf_excluded_fields_do_not_move_identity() {
        let source = vcs_leaf_fixture();
        let base = ledger_vcs_commit_leaf_identity(&source);
        let envelope_a = (1_i64, b"session-a", 10_i64, 20_i64, 1_i64, 1_i64);
        let envelope_b = (2_i64, b"session-b", 30_i64, 40_i64, 2_i64, 2_i64);
        assert_ne!(envelope_a, envelope_b);
        assert_eq!(base, ledger_vcs_commit_leaf_identity(&source));
    }

    #[test]
    fn vcs_commit_root_identity_fields_move_independently() {
        let source = VcsCommitRootIdentitySource {
            leaves: vec![
                identity_test_hash(1),
                identity_test_hash(2),
                identity_test_hash(3),
            ],
        };
        let base = ledger_vcs_commit_root_identity(&source);
        assert_identity_moves(
            "merkle-domain-set",
            &base,
            &ledger_vcs_commit_root_identity_with_domains(
                &source,
                b"frankensim.fs-ledger.vcs.merkle-pair.w2",
                VCS_MERKLE_ODD_PREIMAGE_DOMAIN,
                VCS_COMMIT_ROOT_PREIMAGE_DOMAIN,
            ),
        );
        let mut changed = source.clone();
        changed.leaves.push(identity_test_hash(4));
        for field in ["leaf-count", "tree-shape"] {
            assert_identity_moves(field, &base, &ledger_vcs_commit_root_identity(&changed));
        }
        changed = source.clone();
        changed.leaves.swap(0, 1);
        assert_identity_moves(
            "leaf-order",
            &base,
            &ledger_vcs_commit_root_identity(&changed),
        );
        changed = source;
        changed.leaves[0] = identity_test_hash(5);
        assert_identity_moves(
            "leaf-hashes",
            &base,
            &ledger_vcs_commit_root_identity(&changed),
        );
    }

    #[test]
    fn vcs_commit_root_excluded_fields_do_not_move_identity() {
        let source = VcsCommitRootIdentitySource {
            leaves: vec![identity_test_hash(1), identity_test_hash(2)],
        };
        let base = ledger_vcs_commit_root_identity(&source);
        let envelope_a = (identity_test_hash(3), 1_i64, vec![10_i64, 11], 20_i64);
        let envelope_b = (identity_test_hash(4), 2_i64, vec![30_i64, 31], 40_i64);
        assert_ne!(envelope_a, envelope_b);
        assert_eq!(base, ledger_vcs_commit_root_identity(&source));
    }

    #[test]
    fn vcs_commit_envelope_identity_fields_move_independently() {
        let source = VcsCommitEnvelopeIdentitySource {
            ledger: identity_test_hash(1),
            branch: 2,
            root: identity_test_hash(3),
        };
        let base = ledger_vcs_commit_envelope_identity(&source);
        let mut changed = source.clone();
        changed.ledger = identity_test_hash(4);
        assert_identity_moves(
            "ledger-lineage",
            &base,
            &ledger_vcs_commit_envelope_identity(&changed),
        );
        changed = source.clone();
        changed.branch += 1;
        assert_identity_moves(
            "branch-id",
            &base,
            &ledger_vcs_commit_envelope_identity(&changed),
        );
        changed = source;
        changed.root = identity_test_hash(5);
        assert_identity_moves(
            "semantic-root",
            &base,
            &ledger_vcs_commit_envelope_identity(&changed),
        );
    }

    #[test]
    fn vcs_commit_envelope_excluded_fields_do_not_move_identity() {
        let source = VcsCommitEnvelopeIdentitySource {
            ledger: identity_test_hash(1),
            branch: 2,
            root: identity_test_hash(3),
        };
        let base = ledger_vcs_commit_envelope_identity(&source);
        let envelope_a = (Some(1_i64), Some(identity_test_hash(4)), 5_i64, 6_i64);
        let envelope_b = (Some(7_i64), Some(identity_test_hash(8)), 9_i64, 10_i64);
        assert_ne!(envelope_a, envelope_b);
        assert_eq!(base, ledger_vcs_commit_envelope_identity(&source));
    }

    #[test]
    fn open_migrates_to_current_version() {
        let l = mem();
        assert_eq!(l.schema_version().unwrap(), SCHEMA_VERSION);
        for table in ALL_TABLES {
            // A fresh ledger is empty except the seeded main branch and
            // immutable physical-instance identity rows.
            let expected = u64::from(matches!(*table, "branches" | "ledger_identity"));
            assert_eq!(
                l.table_count(table).unwrap(),
                expected,
                "{table} fresh count"
            );
        }
    }

    #[test]
    fn first_open_rechecks_stale_v0_observation_after_begin() {
        // G4: deterministically place one opener on each side of the
        // preflight-v0 -> BEGIN window. The stale observer must attest and
        // accept the peer's completed initialization, not classify its
        // canonical objects as alien v0 state.
        let path = std::env::temp_dir()
            .join(format!(
                "fs-ledger-first-open-race-{}-{}.db",
                std::process::id(),
                now_wall_ns()
            ))
            .display()
            .to_string();
        let observer_ready = std::sync::Arc::new(std::sync::Barrier::new(2));
        let initializer_done = std::sync::Arc::new(std::sync::Barrier::new(2));
        let observer_path = path.clone();
        let observer_ready_worker = std::sync::Arc::clone(&observer_ready);
        let initializer_done_worker = std::sync::Arc::clone(&initializer_done);
        let observer = std::thread::spawn(move || {
            let conn = Connection::open(&observer_path).expect("open stale observer");
            let ledger = Ledger {
                conn,
                path: observer_path,
                instance_id: LedgerInstanceId([0; 16]),
                read_queries: core::cell::Cell::new(0),
            };
            let observed = ledger.schema_version().expect("preflight user_version");
            assert_eq!(observed, 0, "observer must capture the empty-file marker");
            observer_ready_worker.wait();
            initializer_done_worker.wait();
            ledger
                .migrate_from_observed_version(observed)
                .expect("stale v0 observation accepts canonical peer initialization");
            ledger
                .attest_schema(SCHEMA_VERSION)
                .expect("peer-initialized schema remains exact");
            ledger
                .read_current_instance_id()
                .expect("peer-initialized identity")
        });

        observer_ready.wait();
        let initialized = Ledger::open(&path);
        initializer_done.wait();
        let initialized = initialized.expect("winning first opener");
        let initialized_id = initialized.instance_id();
        assert_eq!(
            observer.join().expect("join stale observer"),
            initialized_id,
            "both first-open contenders must accept one physical authority"
        );
        drop(initialized);
        for suffix in ["", "-wal", "-shm", ".fsqlite-wal", ".fsqlite-shm"] {
            let _ = std::fs::remove_file(format!("{path}{suffix}"));
        }
    }

    #[test]
    fn v8_migration_backfills_verified_claims_and_rejects_corrupt_v7_sources() {
        let (valid, authority) = v7_ledger_with_claim(false);
        valid.migrate().expect("migrate authenticated v7 claim");
        assert_eq!(valid.schema_version().unwrap(), SCHEMA_VERSION);
        assert_eq!(
            valid
                .session_mutation_claim(&authority)
                .expect("read migrated v8 claim")
                .expect("migrated v8 claim")
                .authority,
            authority
        );
        assert_eq!(valid.table_count("session_claim_discovery").unwrap(), 1);

        let (stale_marker, stale_authority) = v7_ledger_with_claim(false);
        for ddl in schema::V8 {
            stale_marker
                .conn
                .execute(ddl)
                .expect("apply exact v8 DDL ahead of stale marker");
        }
        assert_eq!(stale_marker.schema_version().unwrap(), 7);
        stale_marker
            .migrate()
            .expect("heal exact v8 objects with a stale v7 marker");
        assert_eq!(stale_marker.schema_version().unwrap(), SCHEMA_VERSION);
        assert!(
            stale_marker
                .session_mutation_claim(&stale_authority)
                .expect("read stale-marker migrated claim")
                .is_some()
        );

        let (corrupt, _corrupt_authority) = v7_ledger_with_claim(true);
        assert!(matches!(
            corrupt.migrate(),
            Err(LedgerError::Corrupt { .. })
        ));
        assert_eq!(corrupt.schema_version().unwrap(), 7);
        let v8_objects = corrupt
            .conn
            .query(
                "SELECT name FROM sqlite_master \
                 WHERE name = 'session_claim_discovery' LIMIT 1",
            )
            .expect("inspect rolled-back v8 migration");
        assert!(v8_objects.is_empty());
    }

    #[test]
    fn v9_migration_installs_exact_lineage_indexes_and_seals_from_v8() {
        for preapply_v9 in [false, true] {
            let conn = reference_connection(8).expect("construct exact v8 schema");
            // checked_instance_id requires RFC 4122 v4 version/variant bits.
            let mut instance = [0x39; 16];
            instance[6] = 0x49;
            instance[8] = 0x89;
            conn.prepare("INSERT INTO ledger_identity(singleton, instance_id) VALUES (1, ?1)")
                .unwrap()
                .execute_with_params(&[blob_param(&instance)])
                .unwrap();
            if preapply_v9 {
                for ddl in schema::V9 {
                    conn.execute(ddl).expect("preapply exact v9 object");
                }
            }
            conn.execute("PRAGMA user_version = 8").unwrap();
            let ledger = Ledger {
                conn,
                path: ":memory:".to_string(),
                instance_id: LedgerInstanceId(instance),
                read_queries: core::cell::Cell::new(0),
            };
            ledger.migrate().expect("migrate exact v8 to v9");
            assert_eq!(ledger.schema_version().unwrap(), SCHEMA_VERSION);
            assert_eq!(ledger.table_count("artifact_output_seals").unwrap(), 0);
            assert_eq!(ledger.table_count("op_artifact_edge_seals").unwrap(), 0);
            ledger
                .attest_schema(SCHEMA_VERSION)
                .expect("attest v9 schema");
        }
    }

    #[test]
    fn identity_decoder_rejects_hidden_extra_or_noncanonical_rows() {
        assert!(matches!(
            decode_ledger_instance_id(2, None, None),
            Err(LedgerError::InstanceIdentityCorrupt { .. })
        ));
        let uuid = [
            0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x46, 0x77, 0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd,
            0xee, 0xff,
        ];
        let blob = blob_param(&uuid);
        assert!(matches!(
            decode_ledger_instance_id(1, Some(&SqliteValue::Integer(2)), Some(&blob)),
            Err(LedgerError::InstanceIdentityCorrupt { .. })
        ));
    }

    #[test]
    fn signed_database_counts_refuse_negative_values() {
        assert_eq!(nonnegative_u64(0, "fixture").unwrap(), 0);
        assert_eq!(
            nonnegative_u64(i64::MAX, "fixture").unwrap(),
            9_223_372_036_854_775_807u64
        );
        let error = nonnegative_u64(-1, "artifact_info.len").unwrap_err();
        assert!(matches!(
            error,
            LedgerError::Sql { context, detail }
                if context == "artifact_info.len"
                    && detail == "expected non-negative INTEGER, got -1"
        ));
    }

    fn insert_raw_chunked_artifact(
        ledger: &Ledger,
        tag: &[u8],
        len: i64,
        chunk_count: i64,
        chunks: &[(i64, &[u8])],
    ) -> ContentHash {
        let hash = hash_bytes(tag);
        for (seq, bytes) in chunks {
            ledger
                .conn
                .prepare("INSERT INTO artifact_chunks(hash, seq, bytes) VALUES (?1, ?2, ?3)")
                .unwrap()
                .execute_with_params(&[
                    blob_param(hash.as_bytes()),
                    SqliteValue::Integer(*seq),
                    blob_param(bytes),
                ])
                .unwrap();
        }
        ledger
            .conn
            .prepare(
                "INSERT INTO artifacts(hash, kind, bytes, len, chunk_count, meta, created_at) \
                 VALUES (?1, 'corrupt-fixture', NULL, ?2, ?3, NULL, 0)",
            )
            .unwrap()
            .execute_with_params(&[
                blob_param(hash.as_bytes()),
                SqliteValue::Integer(len),
                SqliteValue::Integer(chunk_count),
            ])
            .unwrap();
        hash
    }

    #[test]
    fn inline_artifact_length_mismatch_is_structured_corruption() {
        let ledger = mem();
        let receipt = ledger.put_artifact("blob", b"abc", None).unwrap();
        ledger
            .conn
            .prepare("UPDATE artifacts SET len = 4 WHERE hash = ?1")
            .unwrap()
            .execute_with_params(&[blob_param(receipt.hash.as_bytes())])
            .unwrap();

        let error = ledger.get_artifact(&receipt.hash).unwrap_err();
        assert!(matches!(
            error,
            LedgerError::Corrupt { hash_hex, detail }
                if hash_hex == receipt.hash.to_hex()
                    && detail == "inline length mismatch: recorded 4, found 3"
        ));
    }

    #[test]
    fn chunked_i64_max_length_is_rejected_without_metadata_allocation() {
        let ledger = mem();
        let hash =
            insert_raw_chunked_artifact(&ledger, b"i64-max-length", i64::MAX, 1, &[(0, b"x")]);

        // This metadata used to feed Vec::with_capacity and could panic or
        // attempt an enormous allocation before inspecting the one-byte row.
        let error = ledger.get_artifact(&hash).unwrap_err();
        assert!(matches!(
            error,
            LedgerError::Corrupt { hash_hex, detail }
                if hash_hex == hash.to_hex()
                    && detail
                        == "chunk length mismatch: recorded 9223372036854775807, found 1"
        ));
    }

    #[test]
    fn chunked_count_mismatch_is_rejected_by_metadata_preflight() {
        let ledger = mem();
        let hash = insert_raw_chunked_artifact(&ledger, b"count-mismatch", 1, 2, &[(0, b"x")]);
        let mut prefix = Vec::new();

        let error = ledger
            .read_artifact_chunks(&hash, &mut |chunk| prefix.extend_from_slice(chunk))
            .unwrap_err();
        assert!(prefix.is_empty(), "preflight must run before BLOB delivery");
        assert!(matches!(
            error,
            LedgerError::Corrupt { hash_hex, detail }
                if hash_hex == hash.to_hex()
                    && detail == "chunk count mismatch: recorded 2, found 1"
        ));
    }

    #[test]
    fn chunked_sequence_gap_is_rejected_before_callback() {
        let ledger = mem();
        let hash = insert_raw_chunked_artifact(&ledger, b"sequence-gap", 1, 1, &[(1, b"x")]);
        let mut callback_count = 0;

        let error = ledger
            .read_artifact_chunks(&hash, &mut |_| callback_count += 1)
            .unwrap_err();
        assert_eq!(callback_count, 0, "the invalid row must not be delivered");
        assert!(matches!(
            error,
            LedgerError::Corrupt { hash_hex, detail }
                if hash_hex == hash.to_hex()
                    && detail
                        == "chunk sequence range is not dense from zero: count 1, min 1, max 1"
        ));
    }

    #[test]
    fn chunked_dense_rows_stream_and_materialize_exactly() {
        let ledger = mem();
        let hash = insert_raw_chunked_artifact(&ledger, b"abcde", 5, 2, &[(0, b"abc"), (1, b"de")]);
        let mut chunks = Vec::new();
        let streamed = ledger
            .read_artifact_chunks(&hash, &mut |chunk| chunks.push(chunk.to_vec()))
            .unwrap();

        assert_eq!(streamed, Some(5));
        assert_eq!(chunks, [b"abc".as_slice(), b"de".as_slice()]);
        assert_eq!(ledger.get_artifact(&hash).unwrap().unwrap(), b"abcde");
    }

    #[test]
    fn same_length_inline_tampering_fails_after_streamed_prefix() {
        let ledger = mem();
        let receipt = ledger.put_artifact("blob", b"abc", None).unwrap();
        ledger
            .conn
            .prepare("UPDATE artifacts SET bytes = X'78797A' WHERE hash = ?1")
            .unwrap()
            .execute_with_params(&[blob_param(receipt.hash.as_bytes())])
            .unwrap();

        let mut prefix = Vec::new();
        let error = ledger
            .read_artifact_chunks(&receipt.hash, &mut |bytes| prefix.extend_from_slice(bytes))
            .unwrap_err();
        assert_eq!(
            prefix, b"xyz",
            "late hash failure preserves callback effects"
        );
        assert!(matches!(
            error,
            LedgerError::Corrupt { hash_hex, detail }
                if hash_hex == receipt.hash.to_hex()
                    && detail.starts_with("content hash mismatch: computed ")
        ));
        assert!(matches!(
            ledger.get_artifact(&receipt.hash),
            Err(LedgerError::Corrupt { .. })
        ));
    }

    #[test]
    fn same_length_chunk_tampering_fails_after_streamed_prefix() {
        let ledger = mem();
        let hash = insert_raw_chunked_artifact(&ledger, b"abcde", 5, 2, &[(0, b"abc"), (1, b"de")]);
        ledger
            .conn
            .prepare("UPDATE artifact_chunks SET bytes = X'617863' WHERE hash = ?1 AND seq = 0")
            .unwrap()
            .execute_with_params(&[blob_param(hash.as_bytes())])
            .unwrap();

        let mut prefix = Vec::new();
        let error = ledger
            .read_artifact_chunks(&hash, &mut |bytes| prefix.extend_from_slice(bytes))
            .unwrap_err();
        assert_eq!(
            prefix, b"axcde",
            "late hash failure preserves callback effects"
        );
        assert!(matches!(
            error,
            LedgerError::Corrupt { hash_hex, detail }
                if hash_hex == hash.to_hex()
                    && detail.starts_with("content hash mismatch: computed ")
        ));
    }

    #[test]
    fn storage_preflight_rejects_oversized_blobs_before_callback() {
        let ledger = mem();
        let oversized = vec![0xA5; STORAGE_CHUNK_LEN + 1];
        let inline = ledger.put_artifact("blob", b"inline", None).unwrap();
        ledger
            .conn
            .prepare("UPDATE artifacts SET bytes = ?1, len = ?2 WHERE hash = ?3")
            .unwrap()
            .execute_with_params(&[
                blob_param(&oversized),
                SqliteValue::Integer(int_from_usize(oversized.len())),
                blob_param(inline.hash.as_bytes()),
            ])
            .unwrap();
        let mut callbacks = 0;
        let inline_error = ledger
            .read_artifact_chunks(&inline.hash, &mut |_| callbacks += 1)
            .unwrap_err();
        assert_eq!(callbacks, 0);
        assert!(matches!(
            inline_error,
            LedgerError::Corrupt { detail, .. }
                if detail.starts_with("inline value exceeds the ")
        ));

        let chunked = insert_raw_chunked_artifact(
            &ledger,
            b"oversized-chunk-identity",
            int_from_usize(oversized.len()),
            1,
            &[(0, &oversized)],
        );
        let chunk_error = ledger
            .read_artifact_chunks(&chunked, &mut |_| callbacks += 1)
            .unwrap_err();
        assert_eq!(callbacks, 0);
        assert!(matches!(
            chunk_error,
            LedgerError::Corrupt { detail, .. }
                if detail.starts_with("chunk exceeds the ")
        ));
    }

    #[test]
    fn lint_detects_oversized_and_non_dense_artifact_storage() {
        let ledger = mem();
        let oversized = vec![0x3C; STORAGE_CHUNK_LEN + 1];
        let inline = ledger.put_artifact("blob", b"inline-lint", None).unwrap();
        ledger
            .conn
            .prepare("UPDATE artifacts SET bytes = ?1, len = ?2 WHERE hash = ?3")
            .unwrap()
            .execute_with_params(&[
                blob_param(&oversized),
                SqliteValue::Integer(int_from_usize(oversized.len())),
                blob_param(inline.hash.as_bytes()),
            ])
            .unwrap();
        insert_raw_chunked_artifact(
            &ledger,
            b"oversized-chunk-lint",
            int_from_usize(oversized.len()),
            1,
            &[(0, &oversized)],
        );
        insert_raw_chunked_artifact(&ledger, b"sequence-gap-lint", 1, 1, &[(1, b"x")]);

        let report = ledger.lint().unwrap();
        assert_eq!(report.malformed_artifacts, 2);
        assert_eq!(report.chunk_count_mismatches, 1);
        assert_eq!(report.len_mismatches, 0);
    }

    #[test]
    fn dedupe_refuses_existing_same_length_corruption() {
        let ledger = mem();
        let receipt = ledger.put_artifact("blob", b"abc", None).unwrap();
        ledger
            .conn
            .prepare("UPDATE artifacts SET bytes = X'78797A' WHERE hash = ?1")
            .unwrap()
            .execute_with_params(&[blob_param(receipt.hash.as_bytes())])
            .unwrap();

        assert!(matches!(
            ledger.put_artifact("blob", b"abc", None),
            Err(LedgerError::Corrupt { .. })
        ));
        let mut writer = ledger.artifact_writer("blob").unwrap();
        writer.write(b"abc").unwrap();
        assert!(matches!(
            writer.finish(None),
            Err(LedgerError::Corrupt { .. })
        ));
        assert!(!ledger.in_transaction());
        assert_eq!(ledger.table_count("artifacts").unwrap(), 1);
    }

    #[test]
    fn artifact_dedupe_inline() {
        let l = mem();
        let a = l.put_artifact("blob", b"same bytes", None).unwrap();
        let b = l.put_artifact("blob", b"same bytes", None).unwrap();
        assert!(!a.deduped);
        assert!(b.deduped);
        assert_eq!(a.hash, b.hash);
        assert_eq!(l.table_count("artifacts").unwrap(), 1);
        assert_eq!(l.get_artifact(&a.hash).unwrap().unwrap(), b"same bytes");
    }

    #[test]
    fn bounded_artifact_read_accepts_exact_cap_and_refuses_over_cap_before_callback() {
        let ledger = mem();
        let cap = u64::try_from(MAX_TUNE_MEASURED_BYTES).unwrap();
        let exact_bytes = vec![0xA5; MAX_TUNE_MEASURED_BYTES];
        let exact = ledger
            .put_artifact("bounded-fixture", &exact_bytes, None)
            .unwrap();
        let mut exact_callbacks = 0usize;
        let exact_len = ledger
            .read_artifact_chunks_bounded(&exact.hash, cap, &mut |chunk| {
                exact_callbacks += 1;
                assert_eq!(chunk, exact_bytes.as_slice());
            })
            .unwrap();
        assert_eq!(exact_len, Some(cap));
        assert_eq!(exact_callbacks, 1);
        assert_eq!(
            ledger
                .get_artifact_bounded(&exact.hash, cap)
                .unwrap()
                .as_deref(),
            Some(exact_bytes.as_slice())
        );

        let oversized_bytes = vec![0x5A; MAX_TUNE_MEASURED_BYTES + 1];
        let oversized = ledger
            .put_artifact("bounded-fixture", &oversized_bytes, None)
            .unwrap();
        let mut refused_callbacks = 0usize;
        let error = ledger
            .read_artifact_chunks_bounded(&oversized.hash, cap, &mut |_| {
                refused_callbacks += 1;
            })
            .unwrap_err();
        assert_eq!(refused_callbacks, 0);
        assert!(matches!(
            error,
            LedgerError::ArtifactReadLimit {
                limit,
                observed,
                ..
            } if limit == cap && observed == cap + 1
        ));
        assert!(matches!(
            ledger.get_artifact_bounded(&oversized.hash, cap),
            Err(LedgerError::ArtifactReadLimit {
                limit,
                observed,
                ..
            }) if limit == cap && observed == cap + 1
        ));

        let chunked =
            insert_raw_chunked_artifact(&ledger, b"chunked", 7, 2, &[(0, b"chu"), (1, b"nked")]);
        let mut reconstructed = Vec::new();
        let streamed = ledger
            .read_artifact_chunks_bounded(&chunked, 7, &mut |chunk| {
                reconstructed.extend_from_slice(chunk);
            })
            .unwrap();
        assert_eq!(streamed, Some(7));
        assert_eq!(reconstructed, b"chunked");
        let mut chunked_refused_callbacks = 0usize;
        assert!(matches!(
            ledger.read_artifact_chunks_bounded(&chunked, 6, &mut |_| {
                chunked_refused_callbacks += 1;
            }),
            Err(LedgerError::ArtifactReadLimit {
                limit: 6,
                observed: 7,
                ..
            })
        ));
        assert_eq!(chunked_refused_callbacks, 0);

        let hostile_metadata = ledger
            .put_artifact("bounded-fixture", b"tiny", None)
            .unwrap();
        ledger
            .conn
            .prepare("UPDATE artifacts SET len = ?1 WHERE hash = ?2")
            .unwrap()
            .execute_with_params(&[
                SqliteValue::Integer(int_from_u64(cap + 1)),
                blob_param(hostile_metadata.hash.as_bytes()),
            ])
            .unwrap();
        let mut hostile_callbacks = 0usize;
        assert!(matches!(
            ledger.read_artifact_chunks_bounded(&hostile_metadata.hash, cap, &mut |_| {
                hostile_callbacks += 1;
            }),
            Err(LedgerError::ArtifactReadLimit {
                limit,
                observed,
                ..
            }) if limit == cap && observed == cap + 1
        ));
        assert_eq!(hostile_callbacks, 0);
    }

    #[test]
    fn chunked_duplicate_insert_returns_winner_and_restores_corrupt_prefix() {
        let ledger = mem();
        let hash = insert_raw_chunked_artifact(&ledger, b"abcde", 5, 2, &[(0, b"abc"), (1, b"de")]);
        let outcome = ledger
            .insert_chunked_artifact(&hash, "corrupt-fixture", b"abcde", None)
            .unwrap();
        let ChunkedArtifactInsert::Deduped(info) = outcome else {
            panic!("an existing verified row must be reported as deduped")
        };
        assert_eq!(info.len, 5);
        assert_eq!(info.chunk_count, 2);

        let offered = vec![0x4D; STORAGE_CHUNK_LEN + 1];
        let partial = insert_raw_chunked_artifact(
            &ledger,
            &offered,
            int_from_usize(offered.len()),
            2,
            &[(1, &[0x4D])],
        );
        assert!(matches!(
            ledger.insert_chunked_artifact(&partial, "corrupt-fixture", &offered, None),
            Err(LedgerError::Corrupt { .. })
        ));
        let rows = ledger
            .conn
            .query_with_params(
                "SELECT seq FROM artifact_chunks WHERE hash = ?1 ORDER BY seq",
                &[blob_param(partial.as_bytes())],
            )
            .unwrap();
        assert_eq!(rows.len(), 1, "the attempted seq-0 fill must be removed");
        assert_eq!(row_i64(&rows[0], 0, "restored sequence").unwrap(), 1);
    }

    #[test]
    fn artifact_meta_must_be_json() {
        let l = mem();
        let err = l.put_artifact("blob", b"x", Some("not json")).unwrap_err();
        assert_eq!(err.code(), "LedgerInvalid");
    }

    #[test]
    fn artifact_envelope_write_limits_fail_closed() {
        let ledger = mem();
        let oversized_kind = "k".repeat(MAX_ARTIFACT_KIND_BYTES + 1);
        let error = ledger
            .put_artifact(&oversized_kind, b"kind", None)
            .expect_err("oversized kind must be refused");
        assert!(matches!(
            error,
            LedgerError::Invalid { field, problem }
                if field == "kind" && problem.contains("artifact-kind limit")
        ));

        let oversized_meta = format!("\"{}\"", "m".repeat(MAX_ARTIFACT_META_BYTES));
        let error = ledger
            .put_artifact("blob", b"meta", Some(&oversized_meta))
            .expect_err("oversized metadata must be refused");
        assert!(matches!(
            error,
            LedgerError::Invalid { field, problem }
                if field == "meta" && problem.contains("artifact-metadata limit")
        ));

        let mut writer = ledger.artifact_writer("blob").unwrap();
        writer.write(b"streamed").unwrap();
        assert!(matches!(
            writer.finish(Some(&oversized_meta)),
            Err(LedgerError::Invalid { field, .. }) if field == "meta"
        ));
        assert!(!ledger.in_transaction());
        assert_eq!(ledger.table_count("artifacts").unwrap(), 0);
    }

    #[test]
    fn artifact_info_preflights_hostile_envelope_sizes() {
        let ledger = mem();
        let kind = ledger.put_artifact("blob", b"kind-row", None).unwrap();
        let meta = ledger.put_artifact("blob", b"meta-row", None).unwrap();
        let oversized_kind = "k".repeat(MAX_ARTIFACT_KIND_BYTES + 1);
        let oversized_meta = format!("\"{}\"", "m".repeat(MAX_ARTIFACT_META_BYTES));

        ledger
            .conn
            .prepare("UPDATE artifacts SET kind = ?1 WHERE hash = ?2")
            .unwrap()
            .execute_with_params(&[
                text_param(&oversized_kind),
                blob_param(kind.hash.as_bytes()),
            ])
            .unwrap();
        ledger
            .conn
            .prepare("UPDATE artifacts SET meta = ?1 WHERE hash = ?2")
            .unwrap()
            .execute_with_params(&[
                text_param(&oversized_meta),
                blob_param(meta.hash.as_bytes()),
            ])
            .unwrap();

        assert!(matches!(
            ledger.artifact_info(&kind.hash),
            Err(LedgerError::Corrupt { detail, .. })
                if detail.contains("artifact kind byte length")
        ));
        assert!(matches!(
            ledger.artifact_info(&meta.hash),
            Err(LedgerError::Corrupt { detail, .. })
                if detail.contains("artifact metadata byte length")
        ));
        let lint = ledger.lint().unwrap();
        assert_eq!(lint.malformed_artifacts, 2);
    }

    #[test]
    fn five_explicits_are_enforced_field_by_field() {
        let l = mem();
        let empty_seed = FiveExplicits { seed: &[], ..FX };
        let err = l.begin_op(None, "{}", &empty_seed, 1).unwrap_err();
        assert!(matches!(err, LedgerError::MissingExplicit { ref field, .. } if field == "seed"));

        let bad_budget = FiveExplicits {
            budget: "not json",
            ..FX
        };
        let err = l.begin_op(None, "{}", &bad_budget, 1).unwrap_err();
        assert!(matches!(err, LedgerError::MissingExplicit { ref field, .. } if field == "budget"));

        let empty_versions = FiveExplicits { versions: "", ..FX };
        let err = l.begin_op(None, "{}", &empty_versions, 1).unwrap_err();
        assert!(
            matches!(err, LedgerError::MissingExplicit { ref field, .. } if field == "versions")
        );
    }

    #[test]
    fn op_writes_accept_exact_caps_and_refuse_limit_plus_one_before_json_validation() {
        let ledger = mem();
        let session = vec![0x53; MAX_OP_SESSION_BYTES];
        let seed = vec![0xA7; MAX_OP_SEED_BYTES];
        let ir = json_with_exact_bytes(MAX_OP_IR_BYTES);
        let versions = json_with_exact_bytes(MAX_OP_VERSIONS_BYTES);
        let budget = json_with_exact_bytes(MAX_OP_BUDGET_BYTES);
        let capability = json_with_exact_bytes(MAX_OP_CAPABILITY_BYTES);
        let explicits = FiveExplicits {
            seed: &seed,
            versions: &versions,
            budget: &budget,
            capability: &capability,
        };
        let op = ledger.begin_op(Some(&session), &ir, &explicits, 1).unwrap();
        let diag = json_with_exact_bytes(MAX_OP_DIAG_BYTES);
        ledger.finish_op(op, OpOutcome::Ok, Some(&diag), 2).unwrap();
        let row = ledger.op(op).unwrap().unwrap();
        assert_eq!(row.session.as_deref(), Some(session.as_slice()));
        assert_eq!(row.seed, seed);
        assert_eq!(row.ir.len(), MAX_OP_IR_BYTES);
        assert_eq!(row.versions.len(), MAX_OP_VERSIONS_BYTES);
        assert_eq!(row.budget.len(), MAX_OP_BUDGET_BYTES);
        assert_eq!(row.capability.len(), MAX_OP_CAPABILITY_BYTES);
        assert_eq!(row.diag.as_deref(), Some(diag.as_str()));

        let oversized_session = vec![0; MAX_OP_SESSION_BYTES + 1];
        assert!(matches!(
            ledger.begin_op(Some(&oversized_session), "{}", &FX, 3),
            Err(LedgerError::Invalid { field, problem })
                if field == "session" && problem.contains("operation-field limit")
        ));
        let oversized_ir = "not-json".repeat(MAX_OP_IR_BYTES / 8 + 1);
        assert!(oversized_ir.len() > MAX_OP_IR_BYTES);
        assert!(matches!(
            ledger.begin_op(None, &oversized_ir, &FX, 3),
            Err(LedgerError::Invalid { field, problem })
                if field == "ir" && problem.contains("operation-field limit")
        ));
        let oversized_seed = vec![0; MAX_OP_SEED_BYTES + 1];
        let oversized_seed_fx = FiveExplicits {
            seed: &oversized_seed,
            ..FX
        };
        assert!(matches!(
            ledger.begin_op(None, "{}", &oversized_seed_fx, 3),
            Err(LedgerError::MissingExplicit { field, problem })
                if field == "seed" && problem.contains("operation-field limit")
        ));

        let oversized_json = "not-json".repeat(MAX_OP_FIELD_BYTES / 8 + 1);
        assert!(oversized_json.len() > MAX_OP_FIELD_BYTES);
        for (field, explicits) in [
            (
                "versions",
                FiveExplicits {
                    versions: &oversized_json,
                    ..FX
                },
            ),
            (
                "budget",
                FiveExplicits {
                    budget: &oversized_json,
                    ..FX
                },
            ),
            (
                "capability",
                FiveExplicits {
                    capability: &oversized_json,
                    ..FX
                },
            ),
        ] {
            assert!(matches!(
                ledger.begin_op(None, "{}", &explicits, 3),
                Err(LedgerError::MissingExplicit {
                    field: rejected,
                    problem,
                }) if rejected == field && problem.contains("operation-field limit")
            ));
        }

        let unfinished = ledger.begin_op(None, "{}", &FX, 3).unwrap();
        assert!(matches!(
            ledger.finish_op(unfinished, OpOutcome::Error, Some(&oversized_json), 4),
            Err(LedgerError::Invalid { field, problem })
                if field == "diag" && problem.contains("operation-field limit")
        ));
        assert!(ledger.op(unfinished).unwrap().unwrap().outcome.is_none());
    }

    #[test]
    fn op_reads_preflight_raw_sql_bounds_before_materialization() {
        let ledger = mem();
        let op = ledger.begin_op(None, "{}", &FX, 1).unwrap();
        let exact_ir = json_with_exact_bytes(MAX_OP_IR_BYTES);
        ledger
            .conn
            .prepare("UPDATE ops SET ir = ?1 WHERE id = ?2")
            .unwrap()
            .execute_with_params(&[text_param(&exact_ir), SqliteValue::Integer(op)])
            .unwrap();
        assert_eq!(ledger.op(op).unwrap().unwrap().ir.len(), MAX_OP_IR_BYTES);

        let oversized_ir = json_with_exact_bytes(MAX_OP_IR_BYTES + 1);
        ledger
            .conn
            .prepare("UPDATE ops SET ir = ?1 WHERE id = ?2")
            .unwrap()
            .execute_with_params(&[text_param(&oversized_ir), SqliteValue::Integer(op)])
            .unwrap();
        assert!(matches!(
            ledger.op(op),
            Err(LedgerError::OpCorrupt { op: rejected, .. }) if rejected == op
        ));
        assert_eq!(ledger.lint().unwrap().malformed_ops, 1);

        ledger
            .conn
            .prepare("UPDATE ops SET ir = '{}', versions = ?1 WHERE id = ?2")
            .unwrap()
            .execute_with_params(&[
                text_param(&json_with_exact_bytes(MAX_OP_VERSIONS_BYTES + 1)),
                SqliteValue::Integer(op),
            ])
            .unwrap();
        assert!(matches!(
            ledger.op(op),
            Err(LedgerError::OpCorrupt { op: rejected, .. }) if rejected == op
        ));
        assert_eq!(ledger.lint().unwrap().malformed_ops, 1);

        let oversized_mode = "m".repeat(MAX_OP_FIELD_BYTES + 1);
        ledger
            .conn
            .prepare("UPDATE ops SET versions = '{}', exec_mode = ?1 WHERE id = ?2")
            .unwrap()
            .execute_with_params(&[text_param(&oversized_mode), SqliteValue::Integer(op)])
            .unwrap();
        assert!(matches!(
            ledger.op(op),
            Err(LedgerError::OpCorrupt { op: rejected, .. }) if rejected == op
        ));
        assert_eq!(ledger.lint().unwrap().malformed_ops, 1);
    }

    #[test]
    fn op_execution_context_is_typed_and_rejects_missing_branches() {
        let ledger = mem();
        let main_op = ledger.begin_op(None, "{}", &FX, 1).unwrap();
        assert_eq!(
            ledger.op_execution_context(main_op).unwrap(),
            Some(OpExecutionContext {
                branch: MAIN_BRANCH,
                exec_mode: ExecMode::Deterministic,
            })
        );

        let branch = ledger.fork("fast-context", MAIN_BRANCH).unwrap();
        let fast_op = ledger
            .begin_op_on(branch, ExecMode::Fast, None, "{}", &FX, 2)
            .unwrap();
        assert_eq!(
            ledger.op_execution_context(fast_op).unwrap(),
            Some(OpExecutionContext {
                branch,
                exec_mode: ExecMode::Fast,
            })
        );
        assert_eq!(ledger.op_execution_context(9_999).unwrap(), None);

        ledger
            .conn
            .prepare("UPDATE ops SET branch = 9999 WHERE id = ?1")
            .unwrap()
            .execute_with_params(&[SqliteValue::Integer(fast_op)])
            .unwrap();
        assert!(matches!(
            ledger.op_execution_context(fast_op),
            Err(LedgerError::OpCorrupt { op, detail })
                if op == fast_op && detail.contains("does not exist")
        ));
    }

    #[test]
    fn op_lifecycle_and_double_finish() {
        let l = mem();
        let op = l
            .begin_op(Some(b"sess".as_slice()), r#"{"op":"test"}"#, &FX, 100)
            .unwrap();
        let row = l.op(op).unwrap().unwrap();
        assert_eq!(row.outcome, None);
        assert_eq!(row.seed, FX.seed);
        l.finish_op(op, OpOutcome::Ok, None, 200).unwrap();
        let row = l.op(op).unwrap().unwrap();
        assert_eq!(row.outcome.as_deref(), Some("ok"));
        assert_eq!(row.t_end, Some(200));
        let err = l.finish_op(op, OpOutcome::Error, None, 300).unwrap_err();
        assert_eq!(err, LedgerError::DoubleFinish { op });
        let err = l.finish_op(9999, OpOutcome::Ok, None, 1).unwrap_err();
        assert_eq!(err.code(), "LedgerNotFound");
    }

    #[test]
    fn edges_reject_dangling_references() {
        let l = mem();
        let ghost = hash_bytes(b"never stored");
        let op = l.begin_op(None, "{}", &FX, 1).unwrap();
        let err = l.link(op, &ghost, EdgeRole::Out).unwrap_err();
        assert_eq!(err.code(), "LedgerInvalid");
        let real = l.put_artifact("blob", b"real", None).unwrap();
        l.link(op, &real.hash, EdgeRole::Out).unwrap();
        assert!(l.edge_exists(op, &real.hash, EdgeRole::Out).unwrap());
        assert!(!l.edge_exists(op, &real.hash, EdgeRole::In).unwrap());
        assert!(!l.edge_exists(op + 1, &real.hash, EdgeRole::Out).unwrap());
        assert_eq!(l.table_count("edges").unwrap(), 1);
    }

    #[test]
    fn g0_finish_op_seals_lineage_with_typed_zero_row_refusals() {
        let ledger = mem();
        let retained = ledger.put_artifact("lineage", b"retained", None).unwrap();

        let missing = ledger
            .link(9_999, &retained.hash, EdgeRole::In)
            .unwrap_err();
        assert_eq!(
            missing,
            LedgerError::NotFound {
                what: "op 9999".to_string()
            }
        );

        let corrupt = ledger.begin_op(None, "{}", &FX, 1).unwrap();
        ledger
            .conn
            .prepare("UPDATE ops SET seed = ?1 WHERE id = ?2")
            .unwrap()
            .execute_with_params(&[
                SqliteValue::Blob(vec![0xA5; MAX_OP_SEED_BYTES + 1].into()),
                SqliteValue::Integer(corrupt),
            ])
            .unwrap();
        assert!(matches!(
            ledger.link(corrupt, &retained.hash, EdgeRole::In),
            Err(LedgerError::OpCorrupt { op, .. }) if op == corrupt
        ));
        assert_eq!(ledger.table_count("edges").unwrap(), 0);

        let terminal = ledger.begin_op(None, "{}", &FX, 2).unwrap();
        ledger.link(terminal, &retained.hash, EdgeRole::In).unwrap();
        ledger.finish_op(terminal, OpOutcome::Ok, None, 3).unwrap();
        let late = ledger.put_artifact("lineage", b"late", None).unwrap();
        let edge_count = ledger.table_count("edges").unwrap();

        for (artifact, role) in [(&late.hash, EdgeRole::Out), (&retained.hash, EdgeRole::Out)] {
            let error = ledger.link(terminal, artifact, role).unwrap_err();
            assert_eq!(error, LedgerError::OpLineageSealed { op: terminal });
            assert_eq!(error.code(), "LedgerOpLineageSealed");
        }
        assert_eq!(ledger.table_count("edges").unwrap(), edge_count);
        assert!(
            !ledger
                .edge_exists(terminal, &late.hash, EdgeRole::Out)
                .unwrap()
        );
        assert!(
            !ledger
                .edge_exists(terminal, &retained.hash, EdgeRole::Out)
                .unwrap()
        );
    }

    #[test]
    fn g0_finished_link_refusal_preserves_caller_transaction_rollback() {
        let ledger = mem();
        let retained = ledger
            .put_artifact("transaction-lineage", b"retained", None)
            .unwrap();
        let late = ledger
            .put_artifact("transaction-lineage", b"late", None)
            .unwrap();

        ledger.begin().unwrap();
        let op = ledger.begin_op(None, "{}", &FX, 1).unwrap();
        ledger.link(op, &retained.hash, EdgeRole::Out).unwrap();
        ledger.finish_op(op, OpOutcome::Ok, None, 2).unwrap();
        assert_eq!(
            ledger.link(op, &late.hash, EdgeRole::Out),
            Err(LedgerError::OpLineageSealed { op })
        );
        assert!(ledger.in_transaction());
        ledger.rollback().unwrap();

        assert_eq!(ledger.op(op).unwrap(), None);
        assert!(
            !ledger
                .edge_exists(op, &retained.hash, EdgeRole::Out)
                .unwrap()
        );
        assert!(!ledger.edge_exists(op, &late.hash, EdgeRole::Out).unwrap());
    }

    #[test]
    fn g4_concurrent_link_and_finish_have_one_serial_lineage_order() {
        let path = std::env::temp_dir()
            .join(format!(
                "fs-ledger-link-finish-race-{}-{}.ledger",
                std::process::id(),
                now_wall_ns()
            ))
            .to_string_lossy()
            .into_owned();
        let setup = Ledger::open(&path).unwrap();
        let artifact = setup.put_artifact("race-lineage", b"race", None).unwrap();
        let op = setup.begin_op(None, "{}", &FX, 1).unwrap();
        drop(setup);

        let barrier = std::sync::Arc::new(std::sync::Barrier::new(2));
        let finish_path = path.clone();
        let finish_barrier = barrier.clone();
        let finish = std::thread::spawn(move || {
            let ledger = Ledger::open(&finish_path).unwrap();
            finish_barrier.wait();
            ledger.finish_op(op, OpOutcome::Ok, None, 2)
        });
        let link_path = path.clone();
        let link_hash = artifact.hash;
        let link = std::thread::spawn(move || {
            let ledger = Ledger::open(&link_path).unwrap();
            barrier.wait();
            ledger.link(op, &link_hash, EdgeRole::Out)
        });

        finish.join().unwrap().unwrap();
        let linked_before_finish = match link.join().unwrap() {
            Ok(()) => true,
            Err(LedgerError::OpLineageSealed { op: sealed }) if sealed == op => false,
            Err(error) => panic!("unexpected concurrent link result: {error}"),
        };

        let ledger = Ledger::open(&path).unwrap();
        let row = ledger.op(op).unwrap().expect("terminal raced op");
        assert_eq!(row.outcome.as_deref(), Some("ok"));
        assert_eq!(
            ledger
                .edge_exists(op, &artifact.hash, EdgeRole::Out)
                .unwrap(),
            linked_before_finish,
            "the edge exists iff its atomic insert serialized before finish_op"
        );
        assert_eq!(
            ledger.table_count("edges").unwrap(),
            if linked_before_finish { 1 } else { 0 }
        );
    }

    #[test]
    fn bounded_lineage_queries_report_exact_cap_plus_one_and_roles() {
        let ledger = mem();
        let report = ledger.put_artifact("report", b"report", None).unwrap();
        let other = ledger.put_artifact("other", b"other", None).unwrap();
        let first = ledger.begin_op(None, "{}", &FX, 1).unwrap();
        let second = ledger.begin_op(None, "{}", &FX, 2).unwrap();
        ledger.link(first, &report.hash, EdgeRole::In).unwrap();
        ledger.link(first, &report.hash, EdgeRole::Out).unwrap();
        ledger.link(first, &other.hash, EdgeRole::Out).unwrap();
        ledger.link(second, &report.hash, EdgeRole::Out).unwrap();

        let exact_producers = ledger
            .artifact_producer_ops_bounded(&report.hash, 2)
            .unwrap();
        assert_eq!(exact_producers.op_ids, vec![first, second]);
        assert!(!exact_producers.truncated);
        let capped_producers = ledger
            .artifact_producer_ops_bounded(&report.hash, 1)
            .unwrap();
        assert_eq!(capped_producers.op_ids, vec![first]);
        assert!(capped_producers.truncated);
        let zero_producers = ledger
            .artifact_producer_ops_bounded(&report.hash, 0)
            .unwrap();
        assert!(zero_producers.op_ids.is_empty());
        assert!(zero_producers.truncated);

        let exact_edges = ledger.op_artifact_edges_bounded(first, 3).unwrap();
        assert!(!exact_edges.truncated);
        assert_eq!(exact_edges.edges.len(), 3);
        assert!(exact_edges.edges.contains(&OpArtifactEdge {
            role: EdgeRole::In,
            artifact: report.hash,
        }));
        assert!(exact_edges.edges.contains(&OpArtifactEdge {
            role: EdgeRole::Out,
            artifact: report.hash,
        }));
        assert!(exact_edges.edges.contains(&OpArtifactEdge {
            role: EdgeRole::Out,
            artifact: other.hash,
        }));
        assert_eq!(
            ledger.op_artifact_edges_bounded(first, 3).unwrap(),
            exact_edges,
            "bounded lineage ordering is deterministic"
        );
        let capped_edges = ledger.op_artifact_edges_bounded(first, 2).unwrap();
        assert_eq!(capped_edges.edges.len(), 2);
        assert!(capped_edges.truncated);
        assert_eq!(
            ledger.op_artifact_edges_bounded(9_999, 2).unwrap(),
            BoundedOpArtifactEdges {
                edges: Vec::new(),
                truncated: false,
            }
        );

        for error in [
            ledger
                .artifact_producer_ops_bounded(&report.hash, MAX_LINEAGE_QUERY_ROWS + 1)
                .unwrap_err(),
            ledger
                .op_artifact_edges_bounded(first, MAX_LINEAGE_QUERY_ROWS + 1)
                .unwrap_err(),
        ] {
            assert!(matches!(
                error,
                LedgerError::Invalid { field, problem }
                    if field == "cap" && problem.contains("public maximum")
            ));
        }
    }

    #[test]
    fn lineage_indexes_bound_query_work_and_output_seals_are_immutable() {
        let ledger = mem();
        let report = ledger
            .put_artifact("sealed-report", b"sealed", None)
            .unwrap();
        let source = ledger.put_artifact("source", b"source", None).unwrap();
        let extra = ledger.put_artifact("extra", b"extra", None).unwrap();
        let canonical = ledger.begin_op(None, "{}", &FX, 1).unwrap();
        ledger.link(canonical, &source.hash, EdgeRole::In).unwrap();
        ledger.link(canonical, &report.hash, EdgeRole::Out).unwrap();
        ledger
            .seal_artifact_output(&report.hash, canonical)
            .expect("seal exact sole producer");
        ledger
            .seal_op_artifact_edges(canonical, 2)
            .expect("seal exact two-edge set");
        ledger
            .seal_artifact_output(&report.hash, canonical)
            .expect("same seal is idempotent");
        ledger
            .seal_op_artifact_edges(canonical, 2)
            .expect("same edge-set seal is idempotent");
        assert_eq!(
            ledger.artifact_output_seal(&report.hash).unwrap(),
            Some(canonical)
        );
        assert_eq!(ledger.op_artifact_edge_seal(canonical).unwrap(), Some(2));
        assert!(matches!(
            ledger.link(canonical, &extra.hash, EdgeRole::In),
            Err(LedgerError::Invalid { field, problem })
                if field == "edge" && problem.contains("artifact-edge-set seal")
        ));

        let foreign = ledger.begin_op(None, "{}", &FX, 2).unwrap();
        let error = ledger
            .link(foreign, &report.hash, EdgeRole::Out)
            .expect_err("sealed artifact must reject another producer");
        assert!(matches!(
            error,
            LedgerError::Invalid { field, problem }
                if field == "edge" && problem.contains("exclusive output-producer seal")
        ));
        ledger
            .link(foreign, &report.hash, EdgeRole::In)
            .expect("the output seal does not prohibit input use");
        assert!(matches!(
            ledger.seal_artifact_output(&report.hash, foreign),
            Err(LedgerError::Invalid { field, .. }) if field == "artifact_output_seal"
        ));

        let raw_competing_output = ledger
            .conn
            .prepare("INSERT INTO edges(op, artifact, role) VALUES (?1, ?2, 'out')")
            .unwrap()
            .execute_with_params(&[
                SqliteValue::Integer(foreign),
                blob_param(report.hash.as_bytes()),
            ]);
        assert!(
            raw_competing_output.is_err(),
            "schema trigger is the race backstop"
        );
        for sql in [
            "UPDATE artifact_output_seals SET op = op",
            "DELETE FROM artifact_output_seals",
            "UPDATE op_artifact_edge_seals SET edge_count = edge_count",
            "DELETE FROM op_artifact_edge_seals",
            "DELETE FROM edges WHERE role = 'out'",
        ] {
            assert!(ledger.conn.execute(sql).is_err(), "guard must refuse {sql}");
        }
        assert!(
            ledger
                .conn
                .prepare(
                    "INSERT INTO artifact_output_seals(artifact, op, role) VALUES (?1, ?2, 'out')",
                )
                .unwrap()
                .execute_with_params(&[
                    blob_param(report.hash.as_bytes()),
                    SqliteValue::Integer(canonical),
                ])
                .is_err(),
            "output-seal reinsert guard"
        );
        assert!(
            ledger
                .conn
                .prepare("INSERT INTO op_artifact_edge_seals(op, edge_count) VALUES (?1, 2)",)
                .unwrap()
                .execute_with_params(&[SqliteValue::Integer(canonical)])
                .is_err(),
            "op-edge-seal reinsert guard"
        );

        let producer_plan = ledger
            .conn
            .query_with_params(
                "EXPLAIN QUERY PLAN \
                 SELECT op FROM edges INDEXED BY idx_edges_artifact_role_op \
                 WHERE artifact = ?1 AND role = 'out' ORDER BY op LIMIT 2",
                &[blob_param(report.hash.as_bytes())],
            )
            .unwrap()
            .iter()
            .map(|row| row_text(row, 3, "producer query plan").unwrap())
            .collect::<Vec<_>>()
            .join(" ")
            .to_ascii_lowercase();
        // fsqlite's EXPLAIN QUERY PLAN does not report covering-index usage
        // (it emits "SCAN INDEX ... USING INDEX ...", never sqlite3's
        // "USING COVERING INDEX"), so the covering property itself is not
        // assertable here; the index choice and temp-b-tree absence are.
        assert!(producer_plan.contains("using index idx_edges_artifact_role_op"));
        assert!(!producer_plan.contains("temp b-tree"));

        let edge_plan = ledger
            .conn
            .query_with_params(
                "EXPLAIN QUERY PLAN \
                 SELECT role, artifact FROM edges INDEXED BY idx_edges_op_role_artifact \
                 WHERE op = ?1 ORDER BY role, artifact LIMIT 3",
                &[SqliteValue::Integer(canonical)],
            )
            .unwrap()
            .iter()
            .map(|row| row_text(row, 3, "op-edge query plan").unwrap())
            .collect::<Vec<_>>()
            .join(" ")
            .to_ascii_lowercase();
        assert!(edge_plan.contains("using index idx_edges_op_role_artifact"));
        assert!(!edge_plan.contains("temp b-tree"));
    }

    #[test]
    fn output_seal_read_refuses_a_constraint_bypassed_missing_parent_edge() {
        let ledger = mem();
        let report = ledger
            .put_artifact("sealed-report", b"orphan", None)
            .unwrap();
        let op = ledger.begin_op(None, "{}", &FX, 1).unwrap();
        ledger.link(op, &report.hash, EdgeRole::Out).unwrap();
        ledger.seal_artifact_output(&report.hash, op).unwrap();

        ledger.conn.execute("PRAGMA foreign_keys=OFF").unwrap();
        ledger
            .conn
            .execute("DROP TRIGGER trg_edges_sealed_output_delete")
            .unwrap();
        ledger
            .conn
            .prepare("DELETE FROM edges WHERE op = ?1")
            .unwrap()
            .execute_with_params(&[SqliteValue::Integer(op)])
            .unwrap();
        let delete_guard = schema::V9
            .iter()
            .find(|ddl| ddl.contains("CREATE TRIGGER IF NOT EXISTS trg_edges_sealed_output_delete"))
            .expect("shipped output-edge delete guard");
        ledger.conn.execute(delete_guard).unwrap();
        ledger.conn.execute("PRAGMA foreign_keys=ON").unwrap();

        assert!(matches!(
            ledger.artifact_output_seal(&report.hash),
            Err(LedgerError::Corrupt { hash_hex, detail })
                if hash_hex == report.hash.to_hex() && detail.contains("missing exact output edge")
        ));
        assert!(matches!(
            ledger.seal_artifact_output(&report.hash, op),
            Err(LedgerError::Corrupt { .. })
        ));
    }

    #[test]
    fn op_edge_seal_read_refuses_a_constraint_bypassed_extra_edge() {
        let ledger = mem();
        let first = ledger.put_artifact("edge", b"first", None).unwrap();
        let second = ledger.put_artifact("edge", b"second", None).unwrap();
        let op = ledger.begin_op(None, "{}", &FX, 1).unwrap();
        ledger.link(op, &first.hash, EdgeRole::In).unwrap();
        ledger.seal_op_artifact_edges(op, 1).unwrap();

        ledger
            .conn
            .execute("DROP TRIGGER trg_edges_sealed_op_insert")
            .unwrap();
        ledger
            .conn
            .prepare("INSERT INTO edges(op, artifact, role) VALUES (?1, ?2, 'in')")
            .unwrap()
            .execute_with_params(&[SqliteValue::Integer(op), blob_param(second.hash.as_bytes())])
            .unwrap();
        let insert_guard = schema::V9
            .iter()
            .find(|ddl| ddl.contains("CREATE TRIGGER IF NOT EXISTS trg_edges_sealed_op_insert"))
            .expect("shipped sealed-op insert guard");
        ledger.conn.execute(insert_guard).unwrap();

        assert!(matches!(
            ledger.op_artifact_edge_seal(op),
            Err(LedgerError::OpCorrupt { op: rejected, detail })
                if rejected == op && detail.contains("records 1 edges")
        ));
        assert!(matches!(
            ledger.seal_op_artifact_edges(op, 1),
            Err(LedgerError::OpCorrupt { .. })
        ));
    }

    #[test]
    fn concurrent_output_seal_and_competing_link_cannot_both_commit() {
        let path = std::env::temp_dir()
            .join(format!(
                "fs-ledger-output-seal-race-{}-{}.ledger",
                std::process::id(),
                now_wall_ns()
            ))
            .to_string_lossy()
            .into_owned();
        let setup = Ledger::open(&path).unwrap();
        let report = setup.put_artifact("race-report", b"race", None).unwrap();
        let canonical = setup.begin_op(None, "{}", &FX, 1).unwrap();
        let foreign = setup.begin_op(None, "{}", &FX, 2).unwrap();
        setup.link(canonical, &report.hash, EdgeRole::Out).unwrap();
        drop(setup);

        let barrier = std::sync::Arc::new(std::sync::Barrier::new(2));
        let seal_barrier = barrier.clone();
        let seal_path = path.clone();
        let seal_hash = report.hash;
        let seal = std::thread::spawn(move || {
            let ledger = Ledger::open(&seal_path).unwrap();
            seal_barrier.wait();
            ledger.seal_artifact_output(&seal_hash, canonical).is_ok()
        });
        let link_barrier = barrier;
        let link_path = path.clone();
        let link_hash = report.hash;
        let link = std::thread::spawn(move || {
            let ledger = Ledger::open(&link_path).unwrap();
            link_barrier.wait();
            ledger.link(foreign, &link_hash, EdgeRole::Out).is_ok()
        });
        let sealed = seal.join().unwrap();
        let linked = link.join().unwrap();
        assert_ne!(sealed, linked, "exactly one competing write may commit");

        let ledger = Ledger::open(&path).unwrap();
        let producers = ledger
            .artifact_producer_ops_bounded(&report.hash, 2)
            .unwrap();
        if sealed {
            assert_eq!(producers.op_ids, vec![canonical]);
            assert_eq!(
                ledger.artifact_output_seal(&report.hash).unwrap(),
                Some(canonical)
            );
        } else {
            assert_eq!(producers.op_ids, vec![canonical, foreign]);
            assert_eq!(ledger.artifact_output_seal(&report.hash).unwrap(), None);
        }
    }

    #[test]
    fn concurrent_op_edge_seal_and_third_edge_cannot_both_commit() {
        let path = std::env::temp_dir()
            .join(format!(
                "fs-ledger-op-edge-seal-race-{}-{}.ledger",
                std::process::id(),
                now_wall_ns()
            ))
            .to_string_lossy()
            .into_owned();
        let setup = Ledger::open(&path).unwrap();
        let first = setup.put_artifact("race-edge", b"first", None).unwrap();
        let second = setup.put_artifact("race-edge", b"second", None).unwrap();
        let third = setup.put_artifact("race-edge", b"third", None).unwrap();
        let op = setup.begin_op(None, "{}", &FX, 1).unwrap();
        setup.link(op, &first.hash, EdgeRole::In).unwrap();
        setup.link(op, &second.hash, EdgeRole::Out).unwrap();
        drop(setup);

        let barrier = std::sync::Arc::new(std::sync::Barrier::new(2));
        let seal_barrier = barrier.clone();
        let seal_path = path.clone();
        let seal = std::thread::spawn(move || {
            let ledger = Ledger::open(&seal_path).unwrap();
            seal_barrier.wait();
            ledger.seal_op_artifact_edges(op, 2).is_ok()
        });
        let link_barrier = barrier;
        let link_path = path.clone();
        let third_hash = third.hash;
        let link = std::thread::spawn(move || {
            let ledger = Ledger::open(&link_path).unwrap();
            link_barrier.wait();
            ledger.link(op, &third_hash, EdgeRole::In).is_ok()
        });
        let sealed = seal.join().unwrap();
        let linked = link.join().unwrap();
        assert_ne!(sealed, linked, "exactly one competing write may commit");

        let ledger = Ledger::open(&path).unwrap();
        let edges = ledger.op_artifact_edges_bounded(op, 3).unwrap();
        if sealed {
            assert_eq!(edges.edges.len(), 2);
            assert_eq!(ledger.op_artifact_edge_seal(op).unwrap(), Some(2));
        } else {
            assert_eq!(edges.edges.len(), 3);
            assert_eq!(ledger.op_artifact_edge_seal(op).unwrap(), None);
        }
    }

    #[test]
    fn bounded_op_edges_sanitize_malformed_identities_before_materialization() {
        let ledger = mem();
        let artifact = ledger.put_artifact("edge", b"edge", None).unwrap();
        let op = ledger.begin_op(None, "{}", &FX, 1).unwrap();
        ledger.link(op, &artifact.hash, EdgeRole::In).unwrap();
        ledger.conn.execute("PRAGMA foreign_keys=OFF").unwrap();
        ledger
            .conn
            .prepare("UPDATE edges SET artifact = ?1 WHERE op = ?2")
            .unwrap()
            .execute_with_params(&[
                SqliteValue::Blob(vec![0xA5; MAX_OP_FIELD_BYTES + 1].into()),
                SqliteValue::Integer(op),
            ])
            .unwrap();
        assert!(matches!(
            ledger.op_artifact_edges_bounded(op, 1),
            Err(LedgerError::OpCorrupt { op: rejected, detail })
                if rejected == op && detail.contains("non-blob storage")
        ));
        ledger.conn.execute("PRAGMA foreign_keys=ON").unwrap();
    }

    #[test]
    fn metrics_reject_non_finite() {
        let l = mem();
        let op = l.begin_op(None, "{}", &FX, 1).unwrap();
        assert!(l.record_metric(op, 0, "residual", 1.0e-9).is_ok());
        let err = l.record_metric(op, 1, "residual", f64::NAN).unwrap_err();
        assert_eq!(err.code(), "LedgerInvalid");
    }

    #[test]
    fn events_batch_is_atomic() {
        let l = mem();
        let good = EventRow {
            session: None,
            t: 1,
            kind: "tile_complete",
            payload: None,
        };
        let bad = EventRow {
            session: None,
            t: 2,
            kind: "tile_complete",
            payload: Some("nope"),
        };
        let err = l.append_events(&[good, bad]).unwrap_err();
        assert_eq!(err.code(), "LedgerInvalid");
        assert_eq!(l.table_count("events").unwrap(), 0, "batch rolled back");
        l.append_events(&[good, good]).unwrap();
        assert_eq!(l.table_count("events").unwrap(), 2);
    }

    fn json_with_exact_bytes(len: usize) -> String {
        assert!(len >= 8);
        format!(r#"{{"d":"{}"}}"#, "x".repeat(len - 8))
    }

    fn assert_invalid_field<T>(result: Result<T, LedgerError>, expected_field: &str) {
        let Err(LedgerError::Invalid { field, .. }) = result else {
            panic!("expected LedgerInvalid for {expected_field}");
        };
        assert_eq!(field, expected_field);
    }

    #[test]
    fn tune_upserts() {
        let l = mem();
        l.tune_put(
            "gemm",
            "f64-512",
            b"m1",
            r#"{"mc":256}"#,
            r#"{"gflops":100}"#,
        )
        .unwrap();
        l.tune_put(
            "gemm",
            "f64-512",
            b"m1",
            r#"{"mc":384}"#,
            r#"{"gflops":120}"#,
        )
        .unwrap();
        assert_eq!(l.table_count("tune").unwrap(), 1);
        let row = l.tune_get("gemm", "f64-512", b"m1").unwrap().unwrap();
        assert_eq!(row.params, r#"{"mc":384}"#);
        assert!(l.tune_get("gemm", "f64-512", b"m2").unwrap().is_none());
    }

    #[test]
    fn tune_insert_if_absent_preserves_conflicts_and_transactions() {
        let l = mem();
        let original_params = r#"{"mc":256}"#;
        let original_measured = r#"{"gflops":100}"#;
        l.tune_put_if_absent("gemm", "f64-512", b"m1", original_params, original_measured)
            .expect("insert absent row");
        l.tune_put_if_absent("gemm", "f64-512", b"m1", original_params, original_measured)
            .expect("identical insert is an idempotent no-op");
        l.tune_put_if_absent(
            "gemm",
            "f64-512",
            b"m1",
            r#"{"mc":384}"#,
            r#"{"gflops":120}"#,
        )
        .expect("conflicting insert is a non-overwriting no-op");
        let retained = l
            .tune_get("gemm", "f64-512", b"m1")
            .expect("query")
            .expect("retained row");
        assert_eq!(retained.params, original_params);
        assert_eq!(retained.measured, original_measured);
        assert_eq!(l.table_count("tune").expect("count"), 1);

        let malformed = l
            .tune_put_if_absent("gemm", "bad-json", b"m1", "not-json", "{}")
            .expect_err("malformed params must be refused");
        assert_eq!(malformed.code(), "LedgerInvalid");
        let malformed = l
            .tune_put_if_absent("gemm", "bad-json", b"m1", "{}", "not-json")
            .expect_err("malformed measured evidence must be refused");
        assert_eq!(malformed.code(), "LedgerInvalid");
        assert!(l.tune_get("gemm", "bad-json", b"m1").unwrap().is_none());

        l.begin().expect("begin transaction");
        l.tune_put_if_absent("gemm", "f64-1024", b"m1", "{}", "{}")
            .expect("transactional insert");
        assert!(
            l.tune_get("gemm", "f64-1024", b"m1")
                .expect("query in transaction")
                .is_some()
        );
        l.rollback().expect("rollback transaction");
        assert!(
            l.tune_get("gemm", "f64-1024", b"m1")
                .expect("query after rollback")
                .is_none()
        );
    }

    #[test]
    fn tune_rows_scans_across_machines() {
        let l = mem();
        l.tune_put("axpy", "roofline-v1", b"m1", "{}", r#"{"gbs":50}"#)
            .unwrap();
        l.tune_put("axpy", "roofline-v1", b"m2", "{}", r#"{"gbs":60}"#)
            .unwrap();
        l.tune_put("gemm", "f64-512", b"m1", "{}", "{}").unwrap();
        let rows = l.tune_rows("axpy").unwrap();
        assert_eq!(rows.len(), 2);
        assert!(rows.iter().all(|r| r.kernel == "axpy"));
        assert!(rows.iter().any(|r| r.machine == b"m1"));
        assert!(rows.iter().any(|r| r.machine == b"m2"));
        assert!(l.tune_rows("nonexistent").unwrap().is_empty());
    }

    #[test]
    fn tune_boundaries_accept_exact_limits_and_refuse_limit_plus_one() {
        let l = mem();
        let kernel = "k".repeat(MAX_TUNE_KERNEL_BYTES);
        let shape = "s".repeat(MAX_TUNE_SHAPE_CLASS_BYTES);
        let machine = vec![0_u8; MAX_TUNE_MACHINE_BYTES];
        let params = json_with_exact_bytes(MAX_TUNE_PARAMS_BYTES);
        let measured = json_with_exact_bytes(MAX_TUNE_MEASURED_BYTES);
        l.tune_put(&kernel, &shape, &machine, &params, &measured)
            .expect("every exact tune bound is admitted");
        let stored = l
            .tune_get(&kernel, &shape, &machine)
            .expect("bounded read")
            .expect("exact-limit row");
        assert_eq!(stored.kernel.len(), MAX_TUNE_KERNEL_BYTES);
        assert_eq!(stored.shape_class.len(), MAX_TUNE_SHAPE_CLASS_BYTES);
        assert_eq!(stored.machine, machine);
        assert_eq!(stored.params.len(), MAX_TUNE_PARAMS_BYTES);
        assert_eq!(stored.measured.len(), MAX_TUNE_MEASURED_BYTES);

        let oversized_kernel = "k".repeat(MAX_TUNE_KERNEL_BYTES + 1);
        assert_invalid_field(
            l.tune_put(&oversized_kernel, "s", b"m", "{}", "{}"),
            "kernel",
        );
        let oversized_shape = "s".repeat(MAX_TUNE_SHAPE_CLASS_BYTES + 1);
        assert_invalid_field(
            l.tune_put("k", &oversized_shape, b"m", "{}", "{}"),
            "shape_class",
        );
        let oversized_machine = vec![0_u8; MAX_TUNE_MACHINE_BYTES + 1];
        assert_invalid_field(
            l.tune_put("k", "s", &oversized_machine, "{}", "{}"),
            "machine",
        );
        let oversized_params = json_with_exact_bytes(MAX_TUNE_PARAMS_BYTES + 1);
        assert_invalid_field(
            l.tune_put("k", "s", b"m", &oversized_params, "{}"),
            "params",
        );
        let oversized_measured = json_with_exact_bytes(MAX_TUNE_MEASURED_BYTES + 1);
        assert_invalid_field(
            l.tune_put_if_absent("k", "s", b"m", "{}", &oversized_measured),
            "measured",
        );
    }

    #[test]
    fn tune_identities_refuse_empty_nul_and_noncanonical_bytes() {
        let l = mem();
        for kernel in ["", "contains\0nul", " leading", "unicode-é"] {
            assert_invalid_field(l.tune_put(kernel, "s", b"m", "{}", "{}"), "kernel");
            assert_invalid_field(l.tune_rows(kernel), "kernel");
        }
        for shape in ["", "contains\0nul", "trailing ", "unicode-é"] {
            assert_invalid_field(l.tune_put("k", shape, b"m", "{}", "{}"), "shape_class");
            assert_invalid_field(l.tune_get("k", shape, b"m"), "shape_class");
        }
        assert_invalid_field(l.tune_put("k", "s", b"", "{}", "{}"), "machine");
        assert_invalid_field(l.tune_get("k", "s", b""), "machine");
    }

    #[test]
    fn tune_reads_refuse_oversized_raw_sql_rows_before_payload_materialization() {
        let l = mem();
        let oversized = json_with_exact_bytes(MAX_TUNE_PARAMS_BYTES + 1);
        let insert = l
            .conn
            .prepare(
                "INSERT INTO tune(kernel, shape_class, machine, params, measured) \
                 VALUES (?1, ?2, ?3, ?4, ?5)",
            )
            .unwrap();
        insert
            .execute_with_params(&[
                text_param("raw-kernel"),
                text_param("raw-shape"),
                blob_param(b"raw-machine"),
                text_param(&oversized),
                text_param("{}"),
            ])
            .expect("schema permits valid JSON beyond the API bound");
        insert
            .execute_with_params(&[
                text_param("raw-identity"),
                text_param("bad\0shape"),
                blob_param(b"raw-machine"),
                text_param("{}"),
                text_param("{}"),
            ])
            .expect("schema permits a NUL-bearing shape outside the API contract");
        insert
            .execute_with_params(&[
                text_param("bad\0kernel"),
                text_param("raw-shape"),
                blob_param(b"raw-machine"),
                text_param("{}"),
                text_param("{}"),
            ])
            .expect("schema permits a NUL-bearing kernel outside the API contract");

        assert!(matches!(
            l.tune_get("raw-kernel", "raw-shape", b"raw-machine"),
            Err(LedgerError::TuneCorrupt { detail, .. })
                if detail.contains("params byte length")
        ));
        assert!(matches!(
            l.tune_rows("raw-kernel"),
            Err(LedgerError::TuneCorrupt { detail, .. })
                if detail.contains("params byte length")
        ));
        assert!(matches!(
            l.tune_rows("raw-identity"),
            Err(LedgerError::TuneCorrupt { detail, .. })
                if detail.contains("shape_class is not canonical")
        ));
        assert_eq!(l.lint().unwrap().malformed_tune_rows, 3);
    }

    #[test]
    fn tune_scan_refuses_row_and_aggregate_caps_deterministically() {
        let l = mem();
        let insert = l
            .conn
            .prepare(
                "INSERT INTO tune(kernel, shape_class, machine, params, measured) \
                 VALUES (?1, ?2, ?3, ?4, ?5)",
            )
            .unwrap();
        l.begin().unwrap();
        for index in 0..MAX_TUNE_ROWS_PER_KERNEL {
            insert
                .execute_with_params(&[
                    text_param("row-cap"),
                    text_param(&format!("shape-{index:04}")),
                    blob_param(b"m"),
                    text_param("{}"),
                    text_param("{}"),
                ])
                .unwrap();
        }
        l.commit().unwrap();
        let at_cap = l.tune_rows("row-cap").expect("exact row cap");
        assert_eq!(at_cap.len(), MAX_TUNE_ROWS_PER_KERNEL);
        assert_eq!(at_cap.first().unwrap().shape_class, "shape-0000");
        assert_eq!(
            at_cap.last().unwrap().shape_class,
            format!("shape-{:04}", MAX_TUNE_ROWS_PER_KERNEL - 1)
        );
        insert
            .execute_with_params(&[
                text_param("row-cap"),
                text_param("shape-overflow"),
                blob_param(b"m"),
                text_param("{}"),
                text_param("{}"),
            ])
            .unwrap();
        assert!(matches!(
            l.tune_rows("row-cap"),
            Err(LedgerError::TuneReadLimit {
                resource: "rows",
                limit: MAX_TUNE_ROWS_PER_KERNEL,
                observed_at_least,
                ..
            }) if observed_at_least == MAX_TUNE_ROWS_PER_KERNEL + 1
        ));

        let byte_cap_rows = 256_usize;
        let fixed_row_bytes = "shape-0000".len() + b"m".len() + "{}".len() * 2;
        let kernel_bytes = MAX_TUNE_SCAN_BYTES / byte_cap_rows - fixed_row_bytes;
        let byte_cap_kernel = "k".repeat(kernel_bytes);
        for index in 0..byte_cap_rows {
            insert
                .execute_with_params(&[
                    text_param(&byte_cap_kernel),
                    text_param(&format!("shape-{index:04}")),
                    blob_param(b"m"),
                    text_param("{}"),
                    text_param("{}"),
                ])
                .unwrap();
        }
        let exact_byte_cap = l
            .tune_rows(&byte_cap_kernel)
            .expect("kernel bytes are counted once per exact-cap output row");
        assert_eq!(exact_byte_cap.len(), byte_cap_rows);
        insert
            .execute_with_params(&[
                text_param(&byte_cap_kernel),
                text_param("shape-overflow"),
                blob_param(b"m"),
                text_param("{}"),
                text_param("{}"),
            ])
            .unwrap();
        assert!(matches!(
            l.tune_rows(&byte_cap_kernel),
            Err(LedgerError::TuneReadLimit {
                resource: "materialized_bytes",
                limit: MAX_TUNE_SCAN_BYTES,
                observed_at_least,
                ..
            }) if observed_at_least > MAX_TUNE_SCAN_BYTES
        ));
    }

    #[test]
    fn extension_tables_upsert_and_fetch() {
        let l = mem();
        l.put_extension(
            ExtensionTable::UnsafeCapsules,
            "fs-simd/neon",
            r#"{"lines":134}"#,
        )
        .unwrap();
        l.put_extension(
            ExtensionTable::UnsafeCapsules,
            "fs-simd/neon",
            r#"{"lines":140}"#,
        )
        .unwrap();
        assert_eq!(
            l.get_extension(ExtensionTable::UnsafeCapsules, "fs-simd/neon")
                .unwrap()
                .as_deref(),
            Some(r#"{"lines":140}"#)
        );
        assert!(
            l.get_extension(ExtensionTable::Evidence, "missing")
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn lint_clean_on_fresh_and_populated_ledger() {
        let l = mem();
        assert!(l.lint().unwrap().is_clean());
        let op = l.begin_op(None, "{}", &FX, 1).unwrap();
        let art = l.put_artifact("blob", b"bytes", None).unwrap();
        l.link(op, &art.hash, EdgeRole::Out).unwrap();
        l.record_metric(op, 0, "m", 1.0).unwrap();
        l.finish_op(op, OpOutcome::Ok, None, 2).unwrap();
        assert!(l.lint().unwrap().is_clean());
    }

    #[test]
    fn integrity_detects_corruption() {
        let l = mem();
        let a = l.put_artifact("blob", b"precious bytes", None).unwrap();
        assert!(l.verify_artifact_integrity().unwrap().is_clean());
        l.corrupt_artifact_for_test(&a.hash).unwrap();
        let report = l.verify_artifact_integrity().unwrap();
        assert_eq!(report.corrupted, vec![a.hash.to_hex()]);
    }

    #[test]
    fn writer_streams_and_dedupes() {
        let l = mem();
        let data: Vec<u8> = (0..100_000u32).map(|i| (i % 251) as u8).collect();
        let mut w = l.artifact_writer("field").unwrap();
        for piece in data.chunks(7919) {
            w.write(piece).unwrap();
        }
        let r1 = w.finish(None).unwrap();
        assert_eq!(r1.hash, hash_bytes(&data));
        assert!(!r1.deduped);
        assert_eq!(l.get_artifact(&r1.hash).unwrap().unwrap(), data);
        // Same content again → dedupe, still one artifact row.
        let mut w = l.artifact_writer("field").unwrap();
        w.write(&data).unwrap();
        let r2 = w.finish(None).unwrap();
        assert!(r2.deduped);
        assert_eq!(l.table_count("artifacts").unwrap(), 1);
        assert!(l.lint().unwrap().is_clean());
    }

    #[test]
    fn dropped_writer_leaves_zero_residue() {
        let l = mem();
        let mut w = l.artifact_writer("field").unwrap();
        w.write(&[7u8; 10_000]).unwrap();
        drop(w);
        assert_eq!(l.table_count("artifacts").unwrap(), 0);
        assert_eq!(l.table_count("artifact_chunks").unwrap(), 0);
        assert!(!l.in_transaction());
        assert!(l.lint().unwrap().is_clean());
    }

    #[test]
    fn writer_rejected_inside_transaction() {
        let l = mem();
        l.begin().unwrap();
        let err = l
            .artifact_writer("field")
            .err()
            .map(|e| e.code().to_string());
        assert_eq!(err.as_deref(), Some("LedgerWriterInTransaction"));
        l.rollback().unwrap();
    }

    #[test]
    fn read_query_counter_counts_typed_reads_exactly() {
        // bead vm3i: the counter is the measurable basis for verification
        // query budgets — each typed read API counts exactly once, writes
        // count nothing.
        let l = mem();
        assert_eq!(l.read_queries(), 0);
        l.tune_put("k", "s", b"m", "{}", "{}").unwrap();
        assert_eq!(l.read_queries(), 0, "writes are not read queries");
        let _ = l.tune_get("k", "s", b"m").unwrap();
        assert_eq!(l.read_queries(), 1);
        let _ = l.tune_content_identity("k", "s", b"m").unwrap();
        assert_eq!(l.read_queries(), 2);
        let _ = l.tune_rows("k").unwrap();
        assert_eq!(l.read_queries(), 3);
        let _ = l.op(1).unwrap();
        assert_eq!(l.read_queries(), 4);
        let _ = l.op_execution_context(1).unwrap();
        assert_eq!(l.read_queries(), 5);
        let absent = hash_bytes(b"absent");
        let _ = l.artifact_producer_ops_bounded(&absent, 1).unwrap();
        let _ = l.op_artifact_edges_bounded(1, 1).unwrap();
        assert_eq!(l.read_queries(), 7);
        let _ = l.artifact_output_seal(&absent).unwrap();
        let _ = l.op_artifact_edge_seal(1).unwrap();
        assert_eq!(l.read_queries(), 9);
        let _ = l.get_artifact(&absent).unwrap();
        let _ = l.edge_exists(1, &absent, EdgeRole::Out).unwrap();
        assert_eq!(l.read_queries(), 11);
        let _ = l.checked_instance_id().unwrap();
        assert_eq!(l.read_queries(), 12);
    }
}

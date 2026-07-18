//! Portable semantic bindings for content-addressed runtime-state checkpoints.
//!
//! The generic solver checkpoint proves drain/run provenance. This module
//! instead binds application state to the exact law, law version, state-schema
//! version, canonical parameter block, and injected contract/code identity a
//! replayer says it understands. Stored bytes are returned only after that
//! complete semantic tuple matches.

use std::fmt;
use std::str;

use fs_blake3::hash_domain;
use fsqlite::SqliteValue;

use crate::{ContentHash, Ledger, LedgerError, blob_param, now_wall_ns, sql_err};

/// Semantic version of the portable state-checkpoint receipt identity.
pub const STATE_CHECKPOINT_RECEIPT_IDENTITY_VERSION: u32 = 1;
/// BLAKE3 domain for the portable state-checkpoint receipt identity.
pub const STATE_CHECKPOINT_RECEIPT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-ledger.state-checkpoint-receipt.v1";
/// Exact artifact kind required for canonical constitutive runtime state.
pub const RUNTIME_STATE_ARTIFACT_KIND: &str = "constitutive-runtime-state";
/// Maximum runtime-state bytes materialized by one checkpoint operation.
pub const MAX_RUNTIME_STATE_CHECKPOINT_BYTES: u64 = 64 * 1024 * 1024;
/// Maximum UTF-8 byte length accepted for an fs-matdb law id at the API.
pub const MAX_STATE_CHECKPOINT_LAW_ID_BYTES: usize = 256;

const STATE_CHECKPOINT_RECEIPT_FIXED_TRANSPORT_BYTES: usize =
    4 + 32 + 2 + 4 + 4 + 32 + 32 + 32 + 32;

/// Typed durable identity of one logical state slot.
///
/// The wrapped digest is normally adapted from the Machine-IR `StateSlotId`.
/// Keeping it nominal here prevents implicit exchange with the runtime-state,
/// parameter, or implementation hashes at the persistence boundary without
/// creating an `fs-ledger` -> `fs-ir` dependency cycle. The explicit raw-hash
/// adapter does not itself prove which upstream component minted the digest.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StateSlotId(ContentHash);

impl StateSlotId {
    /// Explicitly adapt a caller-asserted, already domain-separated durable
    /// state-slot digest.
    #[must_use]
    pub const fn from_content_hash(hash: ContentHash) -> Self {
        Self(hash)
    }

    /// Underlying domain-separated state-slot digest.
    #[must_use]
    pub const fn content_hash(self) -> ContentHash {
        self.0
    }
}

/// Exact opaque semantics a caller claims to have available for replay.
///
/// The ledger compares these injected values exactly, but cannot prove which
/// upstream component minted them or interpret an executable implementation.
#[derive(Debug, Clone, Copy)]
pub struct KnownStateSemantics<'a> {
    /// Exact fs-matdb law id.
    pub law_id: &'a str,
    /// Exact law semantic version.
    pub law_version: u32,
    /// Exact runtime-state schema version.
    pub state_schema_version: u32,
    /// Caller-asserted canonical parameter-block identity expected from L1.
    pub canonical_parameters_hash: ContentHash,
    /// Caller-asserted L3 contract plus implementation identity.
    pub contract_and_code_hash: ContentHash,
}

/// Inputs for one immutable semantic state-checkpoint receipt.
#[derive(Debug, Clone, Copy)]
pub struct StateCheckpointClaim<'a> {
    /// Stable logical slot, independent of vector position.
    pub state_slot: StateSlotId,
    /// Complete law/parameter/schema/implementation tuple used to encode the
    /// state. Upstream admission remains responsible for mint provenance.
    pub semantics: KnownStateSemantics<'a>,
    /// Existing `constitutive-runtime-state` artifact.
    pub runtime_state_artifact: ContentHash,
}

/// Portable private-field receipt binding one runtime-state artifact to known
/// semantics.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StateCheckpointReceipt {
    state_slot: StateSlotId,
    law_id: Box<str>,
    law_version: u32,
    state_schema_version: u32,
    runtime_state_artifact: ContentHash,
    canonical_parameters_hash: ContentHash,
    contract_and_code_hash: ContentHash,
    content_hash: ContentHash,
}

/// Exhaustive owner-type classifier for the semantic checkpoint identity.
/// Adding a receipt field must break identity governance until its role is
/// classified deliberately.
#[allow(dead_code)]
fn classify_state_checkpoint_receipt_identity_fields(source: &StateCheckpointReceipt) {
    let StateCheckpointReceipt {
        state_slot,
        law_id,
        law_version,
        state_schema_version,
        runtime_state_artifact,
        canonical_parameters_hash,
        contract_and_code_hash,
        content_hash,
    } = source;
    let _ = (
        state_slot,
        law_id,
        law_version,
        state_schema_version,
        runtime_state_artifact,
        canonical_parameters_hash,
        contract_and_code_hash,
        content_hash,
    );
}

/// Owner-local semantic state-checkpoint declaration consumed by
/// `xtask check-identities`.
#[allow(dead_code)]
pub const STATE_CHECKPOINT_RECEIPT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-ledger:state-checkpoint-receipt",
    "version_const=STATE_CHECKPOINT_RECEIPT_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-ledger.state-checkpoint-receipt.v1",
    "domain_const=STATE_CHECKPOINT_RECEIPT_IDENTITY_DOMAIN",
    "encoder=checkpoint_receipt_hash",
    "encoder_helpers=checkpoint_receipt_hash_with_schema",
    "schema_constants=STATE_CHECKPOINT_RECEIPT_IDENTITY_VERSION,STATE_CHECKPOINT_RECEIPT_IDENTITY_DOMAIN,STATE_CHECKPOINT_RECEIPT_FIXED_TRANSPORT_BYTES,MAX_STATE_CHECKPOINT_LAW_ID_BYTES,RUNTIME_STATE_ARTIFACT_KIND,MAX_RUNTIME_STATE_CHECKPOINT_BYTES,crates/fs-ledger/src/schema.rs#V12",
    "schema_functions=StateCheckpointReceipt::to_bytes,StateCheckpointReceipt::from_bytes,validate_receipt_identity,validate_law_id,receipt_from_semantics,Ledger::stored_state_checkpoint,Ledger::insert_state_checkpoint,Ledger::record_state_checkpoint,Ledger::load_state_checkpoint,Ledger::verify_state_checkpoint_receipt,Ledger::ensure_known_state_semantics,Ledger::load_runtime_state_artifact,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=fs-ledger:artifact-content,fs-matdb:canonical-parameter-block",
    "digest=blake3-256-domain-separated",
    "encoding=typed-binary",
    "sources=StateCheckpointReceipt",
    "source_fields=StateCheckpointReceipt.state_slot:semantic,StateCheckpointReceipt.law_id:semantic,StateCheckpointReceipt.law_version:semantic,StateCheckpointReceipt.state_schema_version:semantic,StateCheckpointReceipt.runtime_state_artifact:semantic,StateCheckpointReceipt.canonical_parameters_hash:semantic,StateCheckpointReceipt.contract_and_code_hash:semantic,StateCheckpointReceipt.content_hash:derived:recomputed-from-semantic-fields",
    "source_bindings=StateCheckpointReceipt.state_slot>state-slot,StateCheckpointReceipt.law_id>law-id-byte-count+law-id-utf8,StateCheckpointReceipt.law_version>law-version,StateCheckpointReceipt.state_schema_version>state-schema-version,StateCheckpointReceipt.runtime_state_artifact>runtime-state-artifact,StateCheckpointReceipt.canonical_parameters_hash>canonical-parameters-hash,StateCheckpointReceipt.contract_and_code_hash>contract-and-code-hash",
    "external_semantic_fields=identity-domain,identity-version,canonical-field-order",
    "semantic_fields=identity-domain,identity-version,canonical-field-order,state-slot,law-id-byte-count,law-id-utf8,law-version,state-schema-version,runtime-state-artifact,canonical-parameters-hash,contract-and-code-hash",
    "excluded_fields=none",
    "consumers=Ledger::record_state_checkpoint,Ledger::load_state_checkpoint,Ledger::verify_state_checkpoint_receipt,StateCheckpointReceipt::from_bytes",
    "mutations=identity-domain:crates/fs-ledger/src/state_checkpoint.rs#receipt_identity_and_transport_bind_every_field,identity-version:crates/fs-ledger/src/state_checkpoint.rs#receipt_identity_and_transport_bind_every_field,canonical-field-order:crates/fs-ledger/src/state_checkpoint.rs#receipt_identity_and_transport_bind_every_field,state-slot:crates/fs-ledger/src/state_checkpoint.rs#receipt_identity_and_transport_bind_every_field,law-id-byte-count:crates/fs-ledger/src/state_checkpoint.rs#receipt_identity_and_transport_bind_every_field,law-id-utf8:crates/fs-ledger/src/state_checkpoint.rs#receipt_identity_and_transport_bind_every_field,law-version:crates/fs-ledger/src/state_checkpoint.rs#receipt_identity_and_transport_bind_every_field,state-schema-version:crates/fs-ledger/src/state_checkpoint.rs#receipt_identity_and_transport_bind_every_field,runtime-state-artifact:crates/fs-ledger/src/state_checkpoint.rs#receipt_identity_and_transport_bind_every_field,canonical-parameters-hash:crates/fs-ledger/src/state_checkpoint.rs#receipt_identity_and_transport_bind_every_field,contract-and-code-hash:crates/fs-ledger/src/state_checkpoint.rs#receipt_identity_and_transport_bind_every_field",
    "nonsemantic_mutations=none",
    "field_guard=classify_state_checkpoint_receipt_identity_fields",
    "transport_guard=validate_receipt_identity",
    "version_guard=crates/fs-ledger/src/state_checkpoint.rs#receipt_identity_and_transport_bind_every_field",
    "coupling_surface=fs-ledger:state-checkpoint-receipt",
];

impl StateCheckpointReceipt {
    /// Stable logical state slot.
    #[must_use]
    pub const fn state_slot(&self) -> StateSlotId {
        self.state_slot
    }

    /// Exact fs-matdb law id.
    #[must_use]
    pub fn law_id(&self) -> &str {
        &self.law_id
    }

    /// Exact law semantic version.
    #[must_use]
    pub const fn law_version(&self) -> u32 {
        self.law_version
    }

    /// Exact runtime-state schema version.
    #[must_use]
    pub const fn state_schema_version(&self) -> u32 {
        self.state_schema_version
    }

    /// Content hash of the retained canonical runtime-state bytes.
    #[must_use]
    pub const fn runtime_state_artifact(&self) -> ContentHash {
        self.runtime_state_artifact
    }

    /// Caller-asserted canonical parameter-block identity.
    #[must_use]
    pub const fn canonical_parameters_hash(&self) -> ContentHash {
        self.canonical_parameters_hash
    }

    /// Caller-asserted L3 contract plus implementation identity.
    #[must_use]
    pub const fn contract_and_code_hash(&self) -> ContentHash {
        self.contract_and_code_hash
    }

    /// Domain-separated receipt identity over every semantic field.
    #[must_use]
    pub const fn content_hash(&self) -> ContentHash {
        self.content_hash
    }

    /// Exact versioned transport for package or process boundaries.
    #[must_use]
    pub fn to_bytes(&self) -> Vec<u8> {
        let law = self.law_id.as_bytes();
        let law_len = u16::try_from(law.len()).unwrap_or(u16::MAX);
        let mut bytes = Vec::with_capacity(
            STATE_CHECKPOINT_RECEIPT_FIXED_TRANSPORT_BYTES.saturating_add(law.len()),
        );
        bytes.extend_from_slice(&STATE_CHECKPOINT_RECEIPT_IDENTITY_VERSION.to_le_bytes());
        bytes.extend_from_slice(self.state_slot.0.as_bytes());
        bytes.extend_from_slice(&law_len.to_le_bytes());
        bytes.extend_from_slice(law);
        bytes.extend_from_slice(&self.law_version.to_le_bytes());
        bytes.extend_from_slice(&self.state_schema_version.to_le_bytes());
        bytes.extend_from_slice(self.runtime_state_artifact.as_bytes());
        bytes.extend_from_slice(self.canonical_parameters_hash.as_bytes());
        bytes.extend_from_slice(self.contract_and_code_hash.as_bytes());
        bytes.extend_from_slice(self.content_hash.as_bytes());
        bytes
    }

    /// Decode a transport candidate without granting ledger membership or
    /// semantic replay authority.
    ///
    /// # Errors
    /// Refuses truncated/extended bytes, future versions, oversized or invalid
    /// law ids, and receipt-hash mismatches.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, StateCheckpointDecodeError> {
        if bytes.len() < STATE_CHECKPOINT_RECEIPT_FIXED_TRANSPORT_BYTES {
            return Err(StateCheckpointDecodeError::Length {
                found: bytes.len(),
                expected: STATE_CHECKPOINT_RECEIPT_FIXED_TRANSPORT_BYTES,
            });
        }
        let version = read_u32(bytes, 0)?;
        if version != STATE_CHECKPOINT_RECEIPT_IDENTITY_VERSION {
            return Err(StateCheckpointDecodeError::UnsupportedVersion { found: version });
        }
        let state_slot = StateSlotId(read_hash(bytes, 4)?);
        let law_len = usize::from(read_u16(bytes, 36)?);
        if law_len == 0 || law_len > MAX_STATE_CHECKPOINT_LAW_ID_BYTES {
            return Err(StateCheckpointDecodeError::LawIdLength {
                found: law_len,
                max: MAX_STATE_CHECKPOINT_LAW_ID_BYTES,
            });
        }
        let expected = STATE_CHECKPOINT_RECEIPT_FIXED_TRANSPORT_BYTES
            .checked_add(law_len)
            .ok_or(StateCheckpointDecodeError::Length {
                found: bytes.len(),
                expected: usize::MAX,
            })?;
        if bytes.len() != expected {
            return Err(StateCheckpointDecodeError::Length {
                found: bytes.len(),
                expected,
            });
        }
        let law_end = 38usize
            .checked_add(law_len)
            .ok_or(StateCheckpointDecodeError::Length {
                found: bytes.len(),
                expected: usize::MAX,
            })?;
        let law_id = str::from_utf8(&bytes[38..law_end])
            .map_err(|_| StateCheckpointDecodeError::LawIdUtf8)?
            .to_string()
            .into_boxed_str();
        let law_version = read_u32(bytes, law_end)?;
        let state_schema_version = read_u32(bytes, law_end + 4)?;
        let runtime_state_artifact = read_hash(bytes, law_end + 8)?;
        let canonical_parameters_hash = read_hash(bytes, law_end + 40)?;
        let contract_and_code_hash = read_hash(bytes, law_end + 72)?;
        let content_hash = read_hash(bytes, law_end + 104)?;
        let receipt = Self {
            state_slot,
            law_id,
            law_version,
            state_schema_version,
            runtime_state_artifact,
            canonical_parameters_hash,
            contract_and_code_hash,
            content_hash,
        };
        validate_receipt_identity(&receipt)?;
        Ok(receipt)
    }
}

/// Bounded transport-decoding refusal for a state-checkpoint receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StateCheckpointDecodeError {
    /// Transport length was truncated, extended, or inconsistent.
    Length { found: usize, expected: usize },
    /// Receipt schema is newer or otherwise unsupported.
    UnsupportedVersion { found: u32 },
    /// Law id length was empty or exceeded the fixed API envelope.
    LawIdLength { found: usize, max: usize },
    /// Law id bytes were not UTF-8.
    LawIdUtf8,
    /// Stored receipt hash did not match the semantic fields.
    HashMismatch,
}

impl fmt::Display for StateCheckpointDecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Length { found, expected } => {
                write!(f, "transport length {found}, expected exactly {expected}")
            }
            Self::UnsupportedVersion { found } => write!(
                f,
                "unsupported state-checkpoint receipt version {found}; expected {}",
                STATE_CHECKPOINT_RECEIPT_IDENTITY_VERSION
            ),
            Self::LawIdLength { found, max } => {
                write!(f, "law id has {found} bytes; require 1..={max}")
            }
            Self::LawIdUtf8 => f.write_str("law id is not UTF-8"),
            Self::HashMismatch => f.write_str("receipt hash does not match semantic fields"),
        }
    }
}

impl std::error::Error for StateCheckpointDecodeError {}

/// State bytes released only after immutable-row, artifact, and caller-known
/// semantic verification.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VerifiedStateCheckpoint {
    receipt: StateCheckpointReceipt,
    state_bytes: Vec<u8>,
}

impl VerifiedStateCheckpoint {
    /// Fully verified portable receipt.
    #[must_use]
    pub const fn receipt(&self) -> &StateCheckpointReceipt {
        &self.receipt
    }

    /// Exact retained runtime-state bytes.
    #[must_use]
    pub fn state_bytes(&self) -> &[u8] {
        &self.state_bytes
    }
}

fn read_u16(bytes: &[u8], offset: usize) -> Result<u16, StateCheckpointDecodeError> {
    let raw = bytes
        .get(offset..offset.saturating_add(2))
        .and_then(|value| value.try_into().ok())
        .ok_or(StateCheckpointDecodeError::Length {
            found: bytes.len(),
            expected: offset.saturating_add(2),
        })?;
    Ok(u16::from_le_bytes(raw))
}

fn read_u32(bytes: &[u8], offset: usize) -> Result<u32, StateCheckpointDecodeError> {
    let raw = bytes
        .get(offset..offset.saturating_add(4))
        .and_then(|value| value.try_into().ok())
        .ok_or(StateCheckpointDecodeError::Length {
            found: bytes.len(),
            expected: offset.saturating_add(4),
        })?;
    Ok(u32::from_le_bytes(raw))
}

fn read_hash(bytes: &[u8], offset: usize) -> Result<ContentHash, StateCheckpointDecodeError> {
    let raw = bytes
        .get(offset..offset.saturating_add(32))
        .and_then(ContentHash::from_slice)
        .ok_or(StateCheckpointDecodeError::Length {
            found: bytes.len(),
            expected: offset.saturating_add(32),
        })?;
    Ok(raw)
}

fn checkpoint_receipt_hash_with_schema(
    receipt: &StateCheckpointReceipt,
    version: u32,
    domain: &str,
) -> ContentHash {
    let law = receipt.law_id.as_bytes();
    let mut preimage = Vec::with_capacity(
        STATE_CHECKPOINT_RECEIPT_FIXED_TRANSPORT_BYTES.saturating_sub(32) + law.len(),
    );
    preimage.extend_from_slice(&version.to_le_bytes());
    preimage.extend_from_slice(receipt.state_slot.0.as_bytes());
    preimage.extend_from_slice(&u16::try_from(law.len()).unwrap_or(u16::MAX).to_le_bytes());
    preimage.extend_from_slice(law);
    preimage.extend_from_slice(&receipt.law_version.to_le_bytes());
    preimage.extend_from_slice(&receipt.state_schema_version.to_le_bytes());
    preimage.extend_from_slice(receipt.runtime_state_artifact.as_bytes());
    preimage.extend_from_slice(receipt.canonical_parameters_hash.as_bytes());
    preimage.extend_from_slice(receipt.contract_and_code_hash.as_bytes());
    hash_domain(domain, &preimage)
}

fn checkpoint_receipt_hash(receipt: &StateCheckpointReceipt) -> ContentHash {
    checkpoint_receipt_hash_with_schema(
        receipt,
        STATE_CHECKPOINT_RECEIPT_IDENTITY_VERSION,
        STATE_CHECKPOINT_RECEIPT_IDENTITY_DOMAIN,
    )
}

fn validate_receipt_identity(
    receipt: &StateCheckpointReceipt,
) -> Result<(), StateCheckpointDecodeError> {
    validate_law_id(receipt.law_id()).map_err(|_| StateCheckpointDecodeError::LawIdLength {
        found: receipt.law_id.len(),
        max: MAX_STATE_CHECKPOINT_LAW_ID_BYTES,
    })?;
    if checkpoint_receipt_hash(receipt) != receipt.content_hash {
        return Err(StateCheckpointDecodeError::HashMismatch);
    }
    Ok(())
}

fn validate_law_id(law_id: &str) -> Result<(), &'static str> {
    if law_id.is_empty() {
        return Err("law id must be nonempty");
    }
    if law_id.len() > MAX_STATE_CHECKPOINT_LAW_ID_BYTES {
        return Err("law id exceeds the 256-byte checkpoint envelope");
    }
    Ok(())
}

fn invalid(field: &str, problem: impl Into<String>) -> LedgerError {
    LedgerError::Invalid {
        field: field.to_string(),
        problem: problem.into(),
    }
}

fn stored_corrupt(slot: StateSlotId, detail: impl Into<String>) -> LedgerError {
    stored_corrupt_hash(slot.0, detail)
}

fn stored_corrupt_hash(identity: ContentHash, detail: impl Into<String>) -> LedgerError {
    LedgerError::Corrupt {
        hash_hex: identity.to_hex(),
        detail: format!(
            "semantic state checkpoint row is corrupt: {}",
            detail.into()
        ),
    }
}

fn hash_from_value(
    value: Option<&SqliteValue>,
    slot: StateSlotId,
    field: &'static str,
) -> Result<ContentHash, LedgerError> {
    let Some(SqliteValue::Blob(bytes)) = value else {
        return Err(stored_corrupt(slot, format!("{field} is not a BLOB")));
    };
    ContentHash::from_slice(bytes).ok_or_else(|| {
        stored_corrupt(
            slot,
            format!(
                "{field} must contain exactly 32 bytes, found {}",
                bytes.len()
            ),
        )
    })
}

fn u32_from_value(
    value: Option<&SqliteValue>,
    slot: StateSlotId,
    field: &'static str,
) -> Result<u32, LedgerError> {
    let Some(SqliteValue::Integer(value)) = value else {
        return Err(stored_corrupt(slot, format!("{field} is not an INTEGER")));
    };
    u32::try_from(*value)
        .map_err(|_| stored_corrupt(slot, format!("{field} is outside the u32 domain")))
}

fn receipt_from_semantics(
    state_slot: StateSlotId,
    semantics: KnownStateSemantics<'_>,
    runtime_state_artifact: ContentHash,
) -> Result<StateCheckpointReceipt, LedgerError> {
    validate_law_id(semantics.law_id)
        .map_err(|problem| invalid("state_checkpoint.law_id", problem))?;
    if state_slot.0 == ContentHash([0; 32]) {
        return Err(invalid(
            "state_checkpoint.state_slot",
            "all-zero state-slot identity is reserved and cannot name durable state",
        ));
    }
    if semantics.canonical_parameters_hash == ContentHash([0; 32]) {
        return Err(invalid(
            "state_checkpoint.canonical_parameters_hash",
            "all-zero parameter-block identity is not fs-matdb authority",
        ));
    }
    if semantics.contract_and_code_hash == ContentHash([0; 32]) {
        return Err(invalid(
            "state_checkpoint.contract_and_code_hash",
            "all-zero implementation identity is unknown semantics",
        ));
    }
    let mut receipt = StateCheckpointReceipt {
        state_slot,
        law_id: semantics.law_id.to_string().into_boxed_str(),
        law_version: semantics.law_version,
        state_schema_version: semantics.state_schema_version,
        runtime_state_artifact,
        canonical_parameters_hash: semantics.canonical_parameters_hash,
        contract_and_code_hash: semantics.contract_and_code_hash,
        content_hash: ContentHash([0; 32]),
    };
    receipt.content_hash = checkpoint_receipt_hash(&receipt);
    Ok(receipt)
}

impl Ledger {
    fn stored_state_checkpoint(
        &self,
        receipt_hash: ContentHash,
    ) -> Result<Option<StateCheckpointReceipt>, LedgerError> {
        const GUARD: &str = "typeof(receipt_hash) = 'blob' AND length(receipt_hash) = 32 \
            AND typeof(state_slot) = 'blob' AND length(state_slot) = 32 \
                AND state_slot != X'0000000000000000000000000000000000000000000000000000000000000000' \
            AND typeof(law_id) = 'blob' AND length(law_id) BETWEEN 1 AND 256 \
            AND typeof(law_version) = 'integer' AND law_version BETWEEN 0 AND 4294967295 \
            AND typeof(state_schema_version) = 'integer' \
                AND state_schema_version BETWEEN 0 AND 4294967295 \
            AND typeof(runtime_state_artifact) = 'blob' \
                AND length(runtime_state_artifact) = 32 \
            AND typeof(canonical_parameters_hash) = 'blob' \
                AND length(canonical_parameters_hash) = 32 \
                AND canonical_parameters_hash != X'0000000000000000000000000000000000000000000000000000000000000000' \
            AND typeof(contract_and_code_hash) = 'blob' \
                AND length(contract_and_code_hash) = 32 \
                AND contract_and_code_hash != X'0000000000000000000000000000000000000000000000000000000000000000' \
            AND typeof(created_at) = 'integer'";
        let query = format!(
            "SELECT \
                CASE WHEN {GUARD} THEN receipt_hash ELSE NULL END, \
                CASE WHEN {GUARD} THEN state_slot ELSE NULL END, \
                CASE WHEN {GUARD} THEN law_id ELSE NULL END, \
                CASE WHEN {GUARD} THEN law_version ELSE NULL END, \
                CASE WHEN {GUARD} THEN state_schema_version ELSE NULL END, \
                CASE WHEN {GUARD} THEN runtime_state_artifact ELSE NULL END, \
                CASE WHEN {GUARD} THEN canonical_parameters_hash ELSE NULL END, \
                CASE WHEN {GUARD} THEN contract_and_code_hash ELSE NULL END \
             FROM semantic_state_checkpoint_receipts \
             WHERE receipt_hash = ?1 LIMIT 2"
        );
        let rows = self
            .conn
            .query_with_params(&query, &[blob_param(receipt_hash.as_bytes())])
            .map_err(|error| sql_err("semantic state checkpoint guarded read", &error))?;
        if rows.is_empty() {
            return Ok(None);
        }
        if rows.len() != 1 {
            return Err(stored_corrupt_hash(
                receipt_hash,
                "one receipt hash names multiple immutable rows",
            ));
        }
        let Some(row) = rows.first() else {
            return Err(stored_corrupt_hash(
                receipt_hash,
                "row disappeared after guarded selection",
            ));
        };
        let content_hash = hash_from_value(row.get(0), StateSlotId(receipt_hash), "receipt_hash")?;
        if content_hash != receipt_hash {
            return Err(stored_corrupt_hash(
                receipt_hash,
                "selected row disagrees with its indexed receipt hash",
            ));
        }
        let state_slot = StateSlotId(hash_from_value(
            row.get(1),
            StateSlotId(receipt_hash),
            "state_slot",
        )?);
        let Some(SqliteValue::Blob(law_id_bytes)) = row.get(2) else {
            return Err(stored_corrupt(state_slot, "law_id is not bounded BLOB"));
        };
        let law_id = str::from_utf8(law_id_bytes)
            .map_err(|_| stored_corrupt(state_slot, "law_id is not canonical UTF-8"))?;
        validate_law_id(law_id)
            .map_err(|problem| stored_corrupt(state_slot, format!("law_id: {problem}")))?;
        let receipt = StateCheckpointReceipt {
            state_slot,
            law_id: law_id.to_string().into_boxed_str(),
            law_version: u32_from_value(row.get(3), state_slot, "law_version")?,
            state_schema_version: u32_from_value(row.get(4), state_slot, "state_schema_version")?,
            runtime_state_artifact: hash_from_value(
                row.get(5),
                state_slot,
                "runtime_state_artifact",
            )?,
            canonical_parameters_hash: hash_from_value(
                row.get(6),
                state_slot,
                "canonical_parameters_hash",
            )?,
            contract_and_code_hash: hash_from_value(
                row.get(7),
                state_slot,
                "contract_and_code_hash",
            )?,
            content_hash,
        };
        validate_receipt_identity(&receipt).map_err(|error| {
            stored_corrupt(
                state_slot,
                format!("receipt identity failed verification: {error}"),
            )
        })?;
        Ok(Some(receipt))
    }

    fn load_runtime_state_artifact(
        &self,
        receipt: &StateCheckpointReceipt,
        stored: bool,
    ) -> Result<Vec<u8>, LedgerError> {
        let refuse = |detail: String| {
            if stored {
                stored_corrupt(receipt.state_slot, detail)
            } else {
                invalid("state_checkpoint.runtime_state_artifact", detail)
            }
        };
        let info = self
            .artifact_info(&receipt.runtime_state_artifact)?
            .ok_or_else(|| {
                refuse(format!(
                    "runtime-state artifact {} does not exist",
                    receipt.runtime_state_artifact
                ))
            })?;
        if info.kind != RUNTIME_STATE_ARTIFACT_KIND {
            return Err(refuse(format!(
                "artifact {} has kind {:?}, require {RUNTIME_STATE_ARTIFACT_KIND:?}",
                receipt.runtime_state_artifact, info.kind
            )));
        }
        self.get_artifact_bounded(
            &receipt.runtime_state_artifact,
            MAX_RUNTIME_STATE_CHECKPOINT_BYTES,
        )?
        .ok_or_else(|| {
            refuse(format!(
                "runtime-state artifact {} disappeared during bounded validation",
                receipt.runtime_state_artifact
            ))
        })
    }

    fn ensure_known_state_semantics(
        &self,
        receipt: &StateCheckpointReceipt,
        known: KnownStateSemantics<'_>,
    ) -> Result<(), LedgerError> {
        let expected =
            receipt_from_semantics(receipt.state_slot, known, receipt.runtime_state_artifact)?;
        let mut differences = Vec::new();
        if receipt.law_id != expected.law_id {
            differences.push("law_id");
        }
        if receipt.law_version != expected.law_version {
            differences.push("law_version");
        }
        if receipt.state_schema_version != expected.state_schema_version {
            differences.push("state_schema_version");
        }
        if receipt.canonical_parameters_hash != expected.canonical_parameters_hash {
            differences.push("canonical_parameters_hash");
        }
        if receipt.contract_and_code_hash != expected.contract_and_code_hash {
            differences.push("contract_and_code_hash");
        }
        if differences.is_empty() {
            return Ok(());
        }
        Err(LedgerError::UnknownStateSemantics {
            state_slot_hex: receipt.state_slot.0.to_hex(),
            stored_law: receipt.law_id.to_string(),
            stored_law_version: receipt.law_version,
            stored_state_schema_version: receipt.state_schema_version,
            expected_law: expected.law_id.to_string(),
            expected_law_version: expected.law_version,
            expected_state_schema_version: expected.state_schema_version,
            stored_parameters_hash: receipt.canonical_parameters_hash.to_hex(),
            expected_parameters_hash: expected.canonical_parameters_hash.to_hex(),
            stored_contract_and_code_hash: receipt.contract_and_code_hash.to_hex(),
            expected_contract_and_code_hash: expected.contract_and_code_hash.to_hex(),
            differences,
        })
    }

    fn insert_state_checkpoint(&self, receipt: &StateCheckpointReceipt) -> Result<(), LedgerError> {
        self.conn
            .prepare(
                "INSERT INTO semantic_state_checkpoint_receipts(\
                    receipt_hash, state_slot, law_id, law_version, state_schema_version, \
                    runtime_state_artifact, canonical_parameters_hash, \
                    contract_and_code_hash, created_at\
                 ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9)",
            )
            .map_err(|error| sql_err("semantic state checkpoint insert prepare", &error))?
            .execute_with_params(&[
                blob_param(receipt.content_hash.as_bytes()),
                blob_param(receipt.state_slot.0.as_bytes()),
                blob_param(receipt.law_id.as_bytes()),
                SqliteValue::Integer(i64::from(receipt.law_version)),
                SqliteValue::Integer(i64::from(receipt.state_schema_version)),
                blob_param(receipt.runtime_state_artifact.as_bytes()),
                blob_param(receipt.canonical_parameters_hash.as_bytes()),
                blob_param(receipt.contract_and_code_hash.as_bytes()),
                SqliteValue::Integer(now_wall_ns()),
            ])
            .map_err(|error| sql_err("semantic state checkpoint insert", &error))?;
        Ok(())
    }

    /// Mint or exactly replay one immutable semantic state-checkpoint receipt.
    ///
    /// The caller injects an opaque canonical parameter hash expected from
    /// fs-matdb and an opaque L3 contract/code identity. The ledger binds those
    /// caller assertions exactly and integrity-checks the retained runtime-state
    /// artifact under a 64 MiB cap before insertion. An exact retry returns the
    /// original receipt. A stable slot may accumulate successive immutable
    /// checkpoints because each distinct state/semantic tuple has a distinct
    /// receipt hash.
    ///
    /// # Errors
    /// Open caller transactions, invalid slot/semantic identities, missing,
    /// wrong-kind, oversized, or corrupt state artifacts, receipt-identity
    /// conflicts, and database failures.
    pub fn record_state_checkpoint(
        &self,
        claim: StateCheckpointClaim<'_>,
    ) -> Result<StateCheckpointReceipt, LedgerError> {
        if self.in_transaction() {
            return Err(invalid(
                "state_checkpoint.transaction",
                "state checkpoint recording must own its transaction; commit or roll back first",
            ));
        }
        let receipt = receipt_from_semantics(
            claim.state_slot,
            claim.semantics,
            claim.runtime_state_artifact,
        )?;
        self.load_runtime_state_artifact(&receipt, false)?;

        self.begin()?;
        let write = (|| match self.stored_state_checkpoint(receipt.content_hash)? {
            Some(stored) if stored == receipt => {
                self.load_runtime_state_artifact(&stored, true)?;
                Ok(stored)
            }
            Some(stored) => Err(invalid(
                "state_checkpoint.content_hash",
                format!(
                    "receipt hash {} resolves to conflicting stored slot {}",
                    receipt.content_hash, stored.state_slot.0
                ),
            )),
            None => {
                self.insert_state_checkpoint(&receipt)?;
                Ok(receipt)
            }
        })();
        match write {
            Ok(receipt) => {
                if let Err(error) = self.commit() {
                    let _ = self.rollback();
                    return Err(error);
                }
                Ok(receipt)
            }
            Err(error) => {
                let _ = self.rollback();
                Err(error)
            }
        }
    }

    /// Load exact runtime-state bytes only under caller-known semantics.
    ///
    /// The receipt hash selects one immutable point in a slot's history.
    /// Absence is `Ok(None)`. A present row is hash-verified, compared against
    /// the exact supplied law/version/state-schema/parameters/code tuple, and
    /// only then is its bounded runtime-state artifact materialized.
    ///
    /// # Errors
    /// Unknown or changed semantics, malformed rows, missing/wrong-kind/
    /// oversized/corrupt artifacts, invalid semantic metadata, and database
    /// failures.
    pub fn load_state_checkpoint(
        &self,
        receipt_hash: ContentHash,
        known: KnownStateSemantics<'_>,
    ) -> Result<Option<VerifiedStateCheckpoint>, LedgerError> {
        let Some(receipt) = self.stored_state_checkpoint(receipt_hash)? else {
            return Ok(None);
        };
        self.ensure_known_state_semantics(&receipt, known)?;
        let state_bytes = self.load_runtime_state_artifact(&receipt, true)?;
        Ok(Some(VerifiedStateCheckpoint {
            receipt,
            state_bytes,
        }))
    }

    /// Re-earn immutable ledger membership, known semantics, and retained
    /// runtime-state integrity for a transport receipt candidate.
    ///
    /// # Errors
    /// Self-inconsistent candidates, missing/conflicting rows, unknown
    /// semantics, invalid semantic metadata, and artifact failures.
    pub fn verify_state_checkpoint_receipt(
        &self,
        receipt: &StateCheckpointReceipt,
        known: KnownStateSemantics<'_>,
    ) -> Result<(), LedgerError> {
        validate_receipt_identity(receipt).map_err(|error| {
            invalid(
                "state_checkpoint.receipt",
                format!("candidate receipt identity refused: {error}"),
            )
        })?;
        let stored = self
            .stored_state_checkpoint(receipt.content_hash)?
            .ok_or_else(|| {
                invalid(
                    "state_checkpoint.content_hash",
                    format!(
                        "receipt {} is not stored by this ledger",
                        receipt.content_hash
                    ),
                )
            })?;
        if stored != *receipt {
            return Err(invalid(
                "state_checkpoint.content_hash",
                "stored state checkpoint differs from the supplied candidate",
            ));
        }
        self.ensure_known_state_semantics(&stored, known)?;
        self.load_runtime_state_artifact(&stored, true).map(|_| ())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{SCHEMA_VERSION, schema};

    fn semantics(code: ContentHash) -> KnownStateSemantics<'static> {
        KnownStateSemantics {
            law_id: "checkpoint-test-law",
            law_version: 3,
            state_schema_version: 7,
            canonical_parameters_hash: ContentHash([0x11; 32]),
            contract_and_code_hash: code,
        }
    }

    fn receipt() -> StateCheckpointReceipt {
        receipt_from_semantics(
            StateSlotId(ContentHash([0x21; 32])),
            semantics(ContentHash([0x41; 32])),
            ContentHash([0x31; 32]),
        )
        .expect("receipt fixture")
    }

    #[test]
    fn receipt_identity_and_transport_bind_every_field() {
        let base = receipt();
        assert_eq!(
            StateCheckpointReceipt::from_bytes(&base.to_bytes()).unwrap(),
            base
        );
        let mut mutations = Vec::new();

        let mut changed = base.clone();
        changed.state_slot = StateSlotId(ContentHash([0x22; 32]));
        mutations.push(changed);
        let mut changed = base.clone();
        changed.law_id = "checkpoint-best-law".into();
        assert_eq!(changed.law_id.len(), base.law_id.len());
        mutations.push(changed);
        let mut changed = base.clone();
        changed.law_id = "short-law".into();
        assert_ne!(changed.law_id.len(), base.law_id.len());
        mutations.push(changed);
        let mut changed = base.clone();
        changed.law_version += 1;
        mutations.push(changed);
        let mut changed = base.clone();
        changed.state_schema_version += 1;
        mutations.push(changed);
        let mut changed = base.clone();
        changed.runtime_state_artifact = ContentHash([0x32; 32]);
        mutations.push(changed);
        let mut changed = base.clone();
        changed.canonical_parameters_hash = ContentHash([0x33; 32]);
        mutations.push(changed);
        let mut changed = base.clone();
        changed.contract_and_code_hash = ContentHash([0x42; 32]);
        mutations.push(changed);

        for changed in mutations {
            assert_ne!(checkpoint_receipt_hash(&changed), base.content_hash);
        }
        assert_ne!(
            checkpoint_receipt_hash_with_schema(
                &base,
                STATE_CHECKPOINT_RECEIPT_IDENTITY_VERSION + 1,
                STATE_CHECKPOINT_RECEIPT_IDENTITY_DOMAIN,
            ),
            base.content_hash
        );
        let foreign_domain_hash = checkpoint_receipt_hash_with_schema(
            &base,
            STATE_CHECKPOINT_RECEIPT_IDENTITY_VERSION,
            "org.frankensim.fs-ledger.state-checkpoint-receipt.foreign",
        );
        assert_ne!(foreign_domain_hash, base.content_hash);
        let mut foreign_domain_receipt = base.clone();
        foreign_domain_receipt.content_hash = foreign_domain_hash;
        assert!(matches!(
            StateCheckpointReceipt::from_bytes(&foreign_domain_receipt.to_bytes()),
            Err(StateCheckpointDecodeError::HashMismatch)
        ));

        let law = base.law_id.as_bytes();
        let mut reordered = Vec::new();
        reordered.extend_from_slice(&STATE_CHECKPOINT_RECEIPT_IDENTITY_VERSION.to_le_bytes());
        reordered.extend_from_slice(base.state_slot.0.as_bytes());
        reordered.extend_from_slice(
            &u16::try_from(law.len())
                .expect("bounded fixture law id")
                .to_le_bytes(),
        );
        reordered.extend_from_slice(law);
        reordered.extend_from_slice(&base.state_schema_version.to_le_bytes());
        reordered.extend_from_slice(&base.law_version.to_le_bytes());
        reordered.extend_from_slice(base.runtime_state_artifact.as_bytes());
        reordered.extend_from_slice(base.canonical_parameters_hash.as_bytes());
        reordered.extend_from_slice(base.contract_and_code_hash.as_bytes());
        let reordered_hash = hash_domain(STATE_CHECKPOINT_RECEIPT_IDENTITY_DOMAIN, &reordered);
        assert_ne!(reordered_hash, base.content_hash);
        let mut reordered_receipt = base.clone();
        reordered_receipt.content_hash = reordered_hash;
        assert!(matches!(
            StateCheckpointReceipt::from_bytes(&reordered_receipt.to_bytes()),
            Err(StateCheckpointDecodeError::HashMismatch)
        ));

        let mut future = base.to_bytes();
        future[..4].copy_from_slice(&2u32.to_le_bytes());
        assert!(matches!(
            StateCheckpointReceipt::from_bytes(&future),
            Err(StateCheckpointDecodeError::UnsupportedVersion { found: 2 })
        ));
        let mut tampered = base.to_bytes();
        tampered[40] ^= 1;
        assert!(matches!(
            StateCheckpointReceipt::from_bytes(&tampered),
            Err(StateCheckpointDecodeError::HashMismatch)
        ));
        let mut malformed_law_length = base.to_bytes();
        let encoded_law_length =
            u16::try_from(base.law_id.len() + 1).expect("bounded fixture law id");
        malformed_law_length[36..38].copy_from_slice(&encoded_law_length.to_le_bytes());
        assert!(matches!(
            StateCheckpointReceipt::from_bytes(&malformed_law_length),
            Err(StateCheckpointDecodeError::Length { .. })
        ));
        let mut invalid_utf8 = base.to_bytes();
        invalid_utf8[38] = 0xff;
        assert!(matches!(
            StateCheckpointReceipt::from_bytes(&invalid_utf8),
            Err(StateCheckpointDecodeError::LawIdUtf8)
        ));
        let mut extended = base.to_bytes();
        extended.push(0);
        assert!(matches!(
            StateCheckpointReceipt::from_bytes(&extended),
            Err(StateCheckpointDecodeError::Length { .. })
        ));
    }

    #[test]
    fn law_id_limit_is_utf8_byte_bounded() {
        let exactly_256_bytes = "é".repeat(128);
        let exact = KnownStateSemantics {
            law_id: &exactly_256_bytes,
            law_version: 3,
            state_schema_version: 7,
            canonical_parameters_hash: ContentHash([0x11; 32]),
            contract_and_code_hash: ContentHash([0x42; 32]),
        };
        let receipt = receipt_from_semantics(
            StateSlotId(ContentHash([0x23; 32])),
            exact,
            ContentHash([0x33; 32]),
        )
        .expect("exact 256-byte UTF-8 law id is admitted");
        assert_eq!(receipt.law_id().len(), MAX_STATE_CHECKPOINT_LAW_ID_BYTES);
        assert_eq!(
            StateCheckpointReceipt::from_bytes(&receipt.to_bytes())
                .expect("exact boundary transport round trips"),
            receipt
        );

        let over_limit = "é".repeat(129);
        let oversized = KnownStateSemantics {
            law_id: &over_limit,
            law_version: exact.law_version,
            state_schema_version: exact.state_schema_version,
            canonical_parameters_hash: exact.canonical_parameters_hash,
            contract_and_code_hash: exact.contract_and_code_hash,
        };
        assert!(matches!(
            receipt_from_semantics(
                StateSlotId(ContentHash([0x24; 32])),
                oversized,
                ContentHash([0x34; 32]),
            ),
            Err(LedgerError::Invalid { field, .. }) if field == "state_checkpoint.law_id"
        ));
    }

    fn drop_v12_objects(ledger: &Ledger) {
        for ddl in [
            "DROP TRIGGER IF EXISTS trg_semantic_state_checkpoint_immutable_receipt_reinsert",
            "DROP TRIGGER IF EXISTS trg_semantic_state_checkpoint_immutable_delete",
            "DROP TRIGGER IF EXISTS trg_semantic_state_checkpoint_immutable_update",
            "DROP INDEX IF EXISTS idx_semantic_state_checkpoint_law",
            "DROP INDEX IF EXISTS idx_semantic_state_checkpoint_slot",
            "DROP INDEX IF EXISTS idx_semantic_state_checkpoint_runtime_artifact",
            "DROP TABLE IF EXISTS semantic_state_checkpoint_receipts",
        ] {
            ledger.conn.execute(ddl).expect("remove v12 fixture object");
        }
    }

    #[test]
    fn genuine_v11_migrates_to_an_empty_v12_table() {
        let ledger = Ledger::open(":memory:").expect("fresh v12 ledger");
        drop_v12_objects(&ledger);
        ledger
            .conn
            .execute("PRAGMA user_version = 11")
            .expect("mark genuine v11 fixture");
        ledger
            .migrate_from_observed_version(11)
            .expect("migrate genuine v11 fixture");
        assert_eq!(ledger.schema_version().unwrap(), SCHEMA_VERSION);
        assert_eq!(
            ledger
                .table_count("semantic_state_checkpoint_receipts")
                .unwrap(),
            0
        );

        ledger
            .conn
            .execute("PRAGMA user_version = 11")
            .expect("install stale v11 marker over exact v12 objects");
        ledger
            .migrate_from_observed_version(11)
            .expect("heal exact pre-applied v12 objects");
        assert_eq!(ledger.schema_version().unwrap(), SCHEMA_VERSION);
    }

    #[test]
    fn divergent_early_v12_object_refuses_before_marker_advances() {
        let ledger = Ledger::open(":memory:").expect("fresh v12 ledger");
        drop_v12_objects(&ledger);
        ledger
            .conn
            .execute("CREATE TABLE semantic_state_checkpoint_receipts(alien INTEGER) STRICT")
            .expect("install divergent early object");
        ledger
            .conn
            .execute("PRAGMA user_version = 11")
            .expect("mark v11 fixture");
        assert!(matches!(
            ledger.migrate_from_observed_version(11),
            Err(LedgerError::SchemaMismatch {
                claimed_version: 11,
                ..
            })
        ));
        assert_eq!(ledger.schema_version().unwrap(), 11);
    }

    #[test]
    fn v12_constraints_guards_and_bounded_reads_fail_closed() {
        let ledger = Ledger::open(":memory:").expect("fresh v12 ledger");
        let runtime = ledger
            .put_artifact(RUNTIME_STATE_ARTIFACT_KIND, b"state", None)
            .expect("runtime state")
            .hash;
        let raw_insert = |receipt: ContentHash,
                          slot: ContentHash,
                          law: &[u8],
                          law_version: i64,
                          state_schema_version: i64,
                          state_artifact: ContentHash,
                          parameters: ContentHash,
                          code: ContentHash| {
            ledger
                .conn
                .prepare(
                    "INSERT INTO semantic_state_checkpoint_receipts(\
                        receipt_hash, state_slot, law_id, law_version, state_schema_version, \
                        runtime_state_artifact, canonical_parameters_hash, \
                        contract_and_code_hash, created_at\
                     ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, 0)",
                )
                .expect("prepare raw v12 row")
                .execute_with_params(&[
                    blob_param(receipt.as_bytes()),
                    blob_param(slot.as_bytes()),
                    blob_param(law),
                    SqliteValue::Integer(law_version),
                    SqliteValue::Integer(state_schema_version),
                    blob_param(state_artifact.as_bytes()),
                    blob_param(parameters.as_bytes()),
                    blob_param(code.as_bytes()),
                ])
        };
        assert!(
            raw_insert(
                ContentHash([0x71; 32]),
                ContentHash([0x72; 32]),
                b"",
                1,
                1,
                runtime,
                ContentHash([0x51; 32]),
                ContentHash([0x61; 32]),
            )
            .is_err()
        );
        assert!(
            raw_insert(
                ContentHash([0x73; 32]),
                ContentHash([0x74; 32]),
                &[b'x'; 257],
                1,
                1,
                runtime,
                ContentHash([0x51; 32]),
                ContentHash([0x61; 32]),
            )
            .is_err()
        );
        assert!(
            raw_insert(
                ContentHash([0x75; 32]),
                ContentHash([0x76; 32]),
                b"law",
                -1,
                1,
                runtime,
                ContentHash([0x51; 32]),
                ContentHash([0x61; 32]),
            )
            .is_err()
        );
        assert!(
            raw_insert(
                ContentHash([0x77; 32]),
                ContentHash([0x78; 32]),
                b"law",
                1,
                4_294_967_296,
                runtime,
                ContentHash([0x51; 32]),
                ContentHash([0x61; 32]),
            )
            .is_err()
        );
        assert!(
            raw_insert(
                ContentHash([0x79; 32]),
                ContentHash([0x7A; 32]),
                b"law",
                1,
                1,
                ContentHash([0x7B; 32]),
                ContentHash([0x51; 32]),
                ContentHash([0x61; 32]),
            )
            .is_err()
        );
        for (receipt, slot, parameters, code) in [
            (
                ContentHash([0x7C; 32]),
                ContentHash([0; 32]),
                ContentHash([0x51; 32]),
                ContentHash([0x61; 32]),
            ),
            (
                ContentHash([0x7D; 32]),
                ContentHash([0x7E; 32]),
                ContentHash([0; 32]),
                ContentHash([0x61; 32]),
            ),
            (
                ContentHash([0x7F; 32]),
                ContentHash([0x80; 32]),
                ContentHash([0x51; 32]),
                ContentHash([0; 32]),
            ),
        ] {
            assert!(
                raw_insert(receipt, slot, b"law", 1, 1, runtime, parameters, code).is_err(),
                "all-zero authority-bearing hashes must fail at storage"
            );
        }
        assert_eq!(
            ledger
                .table_count("semantic_state_checkpoint_receipts")
                .unwrap(),
            0
        );

        let receipt = ledger
            .record_state_checkpoint(StateCheckpointClaim {
                state_slot: StateSlotId(ContentHash([0x81; 32])),
                semantics: semantics(ContentHash([0x82; 32])),
                runtime_state_artifact: runtime,
            })
            .expect("public verified insert");
        let update = ledger
            .conn
            .execute("UPDATE semantic_state_checkpoint_receipts SET created_at = created_at + 1")
            .expect_err("immutable update trigger");
        assert!(
            update
                .to_string()
                .contains("semantic state checkpoint receipt is immutable")
        );
        let delete = ledger
            .conn
            .execute("DELETE FROM semantic_state_checkpoint_receipts")
            .expect_err("immutable delete trigger");
        assert!(
            delete
                .to_string()
                .contains("semantic state checkpoint receipt is immutable")
        );
        let reinsert = ledger
            .conn
            .execute(
                "INSERT INTO semantic_state_checkpoint_receipts(\
                    receipt_hash, state_slot, law_id, law_version, state_schema_version, \
                    runtime_state_artifact, canonical_parameters_hash, \
                    contract_and_code_hash, created_at\
                 ) SELECT receipt_hash, state_slot, law_id, law_version, \
                          state_schema_version, runtime_state_artifact, \
                          canonical_parameters_hash, contract_and_code_hash, created_at \
                   FROM semantic_state_checkpoint_receipts LIMIT 1",
            )
            .expect_err("immutable receipt reinsert trigger");
        assert!(
            reinsert
                .to_string()
                .contains("semantic state checkpoint receipt is immutable")
        );
        assert_eq!(
            ledger
                .load_state_checkpoint(receipt.content_hash(), semantics(ContentHash([0x82; 32])),)
                .unwrap()
                .unwrap()
                .state_bytes(),
            b"state"
        );

        let hostile = ContentHash([0x91; 32]);
        raw_insert(
            hostile,
            ContentHash([0x92; 32]),
            &[0xFF],
            1,
            1,
            runtime,
            ContentHash([0x51; 32]),
            ContentHash([0x61; 32]),
        )
        .expect("inject bounded non-UTF8 raw row");
        assert!(matches!(
            ledger.load_state_checkpoint(
                hostile,
                semantics(ContentHash([0x82; 32])),
            ),
            Err(LedgerError::Corrupt { detail, .. }) if detail.contains("UTF-8")
        ));
    }

    #[test]
    fn oversized_runtime_state_refuses_before_receipt_publication() {
        let ledger = Ledger::open(":memory:").expect("fresh v12 ledger");
        let state = ledger
            .put_artifact(RUNTIME_STATE_ARTIFACT_KIND, b"tiny", None)
            .expect("tiny state artifact")
            .hash;
        ledger
            .conn
            .prepare("UPDATE artifacts SET len = ?1 WHERE hash = ?2")
            .expect("prepare hostile metadata fixture")
            .execute_with_params(&[
                SqliteValue::Integer(
                    i64::try_from(MAX_RUNTIME_STATE_CHECKPOINT_BYTES + 1).unwrap(),
                ),
                blob_param(state.as_bytes()),
            ])
            .expect("inflate declared state length");
        assert!(matches!(
            ledger.record_state_checkpoint(StateCheckpointClaim {
                state_slot: StateSlotId(ContentHash([0xA1; 32])),
                semantics: semantics(ContentHash([0xA2; 32])),
                runtime_state_artifact: state,
            }),
            Err(LedgerError::ArtifactReadLimit {
                limit: MAX_RUNTIME_STATE_CHECKPOINT_BYTES,
                observed,
                ..
            }) if observed == MAX_RUNTIME_STATE_CHECKPOINT_BYTES + 1
        ));
        assert_eq!(
            ledger
                .table_count("semantic_state_checkpoint_receipts")
                .unwrap(),
            0
        );
    }

    #[test]
    fn migration_ladder_preserves_v12_through_v16_before_the_v17_batch() {
        assert_eq!(
            schema::MIGRATIONS.len(),
            usize::try_from(SCHEMA_VERSION).unwrap()
        );
        assert_eq!(schema::MIGRATIONS.get(11), Some(&schema::V12));
        assert_eq!(schema::MIGRATIONS.get(12), Some(&schema::V13));
        assert_eq!(schema::MIGRATIONS.get(13), Some(&schema::V14));
        assert_eq!(schema::MIGRATIONS.get(14), Some(&schema::V15));
        assert_eq!(schema::MIGRATIONS.get(15), Some(&schema::V16));
        assert_eq!(schema::MIGRATIONS.last(), Some(&schema::V17));
    }
}

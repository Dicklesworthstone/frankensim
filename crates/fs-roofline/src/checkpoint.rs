//! Bounded, authenticated staleness-history checkpoints (bead vm3i).
//!
//! The exhaustive verifier ([`crate::staleness_at`]) re-validates every
//! matching tune row on every probe: op fetch, protocol parse, dependency
//! receipt, manifest membership, then a per-sibling `tune_get` walk — for H
//! retained runs of a K-kernel registry that is O(H·(K+4)) ledger queries,
//! and sustained CI turns it into an unbounded resource sink.
//!
//! A [`StalenessCheckpoint`] is one exhaustive pass SEALED: per covered row
//! it retains a content hash (over key + params + measured bytes), the
//! build identity, the op-bound dependency-receipt digests, the recorded
//! timestamp, and the verdict. Checkpoints live in the ledger's own tune
//! table under a reserved kernel name, machine-keyed exactly like the
//! production rows they cover, and chain through a domain-separated digest:
//! `digest_i = H(domain, prev_digest_bytes ‖ body_bytes)`. Authentication
//! is the chain over verifier-produced content — operator-observed,
//! tamper-EVIDENT, no cryptographic unforgeability claim (the workspace
//! no-crypto no-claim applies: anyone who can write tune rows can mint a
//! chain; what they cannot do is alter covered history under an EXISTING
//! chain without detection).
//!
//! The fast path ([`staleness_at_checkpointed`]) costs a constant number of
//! ledger queries plus the exhaustive cost of only the DELTA (rows newer
//! than the checkpoint):
//! 1. fetch checkpoint rows (1 bounded query), verify ordinal contiguity
//!    and the digest chain — any break FAILS CLOSED to the exhaustive path;
//! 2. fetch the kernel's rows (1 bounded query) and classify the cheap
//!    lattice prefix identically to the exhaustive path;
//! 3. every covered row must hash to its checkpointed entry (tamper and
//!    rollback → `CorruptEvidence`); every checkpointed row must still
//!    exist (removal → `CorruptEvidence`); corrupt verdicts are TOMBSTONES
//!    (once corrupt, forever corrupt — re-checkpointing preserves them, so
//!    compaction can refuse but never un-corrupt);
//! 4. covered valid rows replay the build/dependency scan from checkpointed
//!    metadata against the CURRENT build and dependency binding — no op
//!    fetches; delta rows run the full exhaustive validator.

use fs_ledger::{Ledger, LedgerError};

use crate::{
    BuildRowScan, DependencyReceiptBinding, RowSelection, Staleness, ValidatedRooflineRow,
    classify_scanned_rows, executable_build_identity, roofline_machine_key, select_matching_rows,
    validate_roofline_row,
};

/// Reserved tune-table kernel prefix for checkpoint rows. Distinct from any
/// production kernel name, so production row queries never see checkpoints
/// and vice versa.
const CHECKPOINT_KERNEL_PREFIX: &str = "roofline-staleness-checkpoint:";
const CHECKPOINT_SHAPE_PREFIX: &str = "roofline-ckpt-v2:";
const CHECKPOINT_SCHEMA: &str = "fs-roofline-staleness-checkpoint-v2";
const LEGACY_CHECKPOINT_SHAPE_PREFIX: &str = "roofline-ckpt-v1:";
const LEGACY_CHECKPOINT_SCHEMA: &str = "fs-roofline-staleness-checkpoint-v1";
const LEGACY_CHECKPOINT_CHAIN_DOMAIN: &str =
    "org.frankensim.fs-roofline.staleness-checkpoint-chain.v1";
/// Semantic version of the durable staleness checkpoint-chain receipt.
pub const STALENESS_CHECKPOINT_CHAIN_IDENTITY_VERSION: u32 = 2;
/// BLAKE3 derive-key domain for the durable staleness checkpoint-chain receipt.
pub const CHECKPOINT_CHAIN_DOMAIN: &str =
    "org.frankensim.fs-roofline.staleness-checkpoint-chain.v2";
/// Semantic version of the durable tune-row content receipt.
pub const STALENESS_ROW_CONTENT_IDENTITY_VERSION: u32 = 1;
/// BLAKE3 derive-key domain for the durable tune-row content receipt.
pub const ROW_CONTENT_DOMAIN: &str = "org.frankensim.fs-roofline.staleness-row-content.v1";

/// Owner-local tune-row content declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const STALENESS_ROW_CONTENT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-roofline:staleness-row-content",
    "version_const=STALENESS_ROW_CONTENT_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-roofline.staleness-row-content.v1",
    "domain_const=ROW_CONTENT_DOMAIN",
    "encoder=staleness_row_content_receipt",
    "encoder_helpers=StalenessRowContentIdentityInput::from_row,staleness_row_content_receipt_with_domain,push_row_content_field,row_content_hash",
    "schema_constants=STALENESS_ROW_CONTENT_IDENTITY_VERSION,ROW_CONTENT_DOMAIN",
    "schema_functions=seal_entry,row_content_is_current,staleness_at_checkpointed_with,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=none",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=StalenessRowContentIdentityInput",
    "source_fields=StalenessRowContentIdentityInput.kernel:semantic,StalenessRowContentIdentityInput.shape_class:semantic,StalenessRowContentIdentityInput.machine:semantic,StalenessRowContentIdentityInput.params:semantic,StalenessRowContentIdentityInput.measured:semantic",
    "source_bindings=StalenessRowContentIdentityInput.kernel>kernel-bytes,StalenessRowContentIdentityInput.shape_class>shape-class-bytes,StalenessRowContentIdentityInput.machine>machine-key-bytes,StalenessRowContentIdentityInput.params>params-bytes,StalenessRowContentIdentityInput.measured>measured-bytes",
    "external_semantic_fields=digest-domain,identity-version,u64-length-prefix-framing",
    "semantic_fields=digest-domain,identity-version,u64-length-prefix-framing,kernel-bytes,shape-class-bytes,machine-key-bytes,params-bytes,measured-bytes",
    "excluded_fields=none",
    "consumers=CheckpointEntry::row_hash,seal_entry,row_content_is_current,staleness_at_checkpointed_with",
    "mutations=digest-domain:crates/fs-roofline/src/checkpoint.rs#staleness_row_content_identity_fields_move_independently,identity-version:crates/fs-roofline/src/checkpoint.rs#staleness_row_content_identity_versions_fail_closed,u64-length-prefix-framing:crates/fs-roofline/src/checkpoint.rs#staleness_row_content_identity_fields_move_independently,kernel-bytes:crates/fs-roofline/src/checkpoint.rs#staleness_row_content_identity_fields_move_independently,shape-class-bytes:crates/fs-roofline/src/checkpoint.rs#staleness_row_content_identity_fields_move_independently,machine-key-bytes:crates/fs-roofline/src/checkpoint.rs#staleness_row_content_identity_fields_move_independently,params-bytes:crates/fs-roofline/src/checkpoint.rs#staleness_row_content_identity_fields_move_independently,measured-bytes:crates/fs-roofline/src/checkpoint.rs#staleness_row_content_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_staleness_row_content_identity_fields",
    "transport_guard=row_content_is_current",
    "version_guard=crates/fs-roofline/src/checkpoint.rs#staleness_row_content_identity_versions_fail_closed",
    "coupling_surface=fs-roofline:staleness-row-content",
];

/// Owner-local checkpoint-chain declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const STALENESS_CHECKPOINT_CHAIN_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-roofline:staleness-checkpoint-chain",
    "version_const=STALENESS_CHECKPOINT_CHAIN_IDENTITY_VERSION",
    "version=2",
    "domain=org.frankensim.fs-roofline.staleness-checkpoint-chain.v2",
    "domain_const=CHECKPOINT_CHAIN_DOMAIN",
    "encoder=staleness_checkpoint_chain_receipt",
    "encoder_helpers=staleness_checkpoint_chain_receipt_with_domain,checkpoint_chain_body_json,checkpoint_chain_body_json_with_schema,entry_json,lower_hex,decode_lower_hex_utf8",
    "schema_constants=STALENESS_CHECKPOINT_CHAIN_IDENTITY_VERSION,CHECKPOINT_CHAIN_DOMAIN,CHECKPOINT_SCHEMA,CHECKPOINT_KERNEL_PREFIX,CHECKPOINT_SHAPE_PREFIX,LEGACY_CHECKPOINT_SHAPE_PREFIX,LEGACY_CHECKPOINT_SCHEMA,LEGACY_CHECKPOINT_CHAIN_DOMAIN",
    "schema_functions=checkpoint_kernel,checkpoint_shape_class,checkpoint_params,parse_body,legacy_checkpoint_shape_class,legacy_checkpoint_params,legacy_entry_json,legacy_body_json,legacy_chain_digest,parse_legacy_body,normalize_legacy_entries,verify_legacy_chain,verify_current_chain,load_chain,checkpoint_staleness_history,staleness_at_checkpointed_with,crates/fs-blake3/src/lib.rs#ContentHash::as_bytes,crates/fs-blake3/src/lib.rs#ContentHash::from_hex,crates/fs-blake3/src/lib.rs#ContentHash::to_hex,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=fs-roofline:staleness-row-content,fs-roofline:executable-build,fs-la:depgraph-receipt",
    "digest=fs-blake3",
    "encoding=canonical-transport-exact-bits",
    "sources=StalenessCheckpointChainIdentityInput,CheckpointEntry",
    "source_fields=StalenessCheckpointChainIdentityInput.kernel:semantic,StalenessCheckpointChainIdentityInput.version:semantic,StalenessCheckpointChainIdentityInput.ordinal:semantic,StalenessCheckpointChainIdentityInput.previous:semantic,StalenessCheckpointChainIdentityInput.entries:semantic,CheckpointEntry.shape_class:semantic,CheckpointEntry.row_hash:semantic,CheckpointEntry.build:semantic,CheckpointEntry.dep_digest:semantic,CheckpointEntry.dep_artifact:semantic,CheckpointEntry.recorded_at_ns:semantic,CheckpointEntry.verdict:semantic",
    "source_bindings=StalenessCheckpointChainIdentityInput.kernel>kernel,StalenessCheckpointChainIdentityInput.version>kernel-version,StalenessCheckpointChainIdentityInput.ordinal>checkpoint-ordinal,StalenessCheckpointChainIdentityInput.previous>previous-digest-presence+previous-digest-value,StalenessCheckpointChainIdentityInput.entries>entry-count+ordered-entries,CheckpointEntry.shape_class>entry-shape-class,CheckpointEntry.row_hash>row-content-child,CheckpointEntry.build>executable-build-child,CheckpointEntry.dep_digest>dependency-domain-digest,CheckpointEntry.dep_artifact>dependency-artifact-hash,CheckpointEntry.recorded_at_ns>recorded-wall-ns,CheckpointEntry.verdict>verdict-tag",
    "external_semantic_fields=digest-domain,identity-version,checkpoint-body-schema,canonical-json-layout,lowercase-hex-string-framing,previous-digest-framing,legacy-v1-chain-binding",
    "semantic_fields=digest-domain,identity-version,checkpoint-body-schema,canonical-json-layout,lowercase-hex-string-framing,previous-digest-framing,legacy-v1-chain-binding,kernel,kernel-version,checkpoint-ordinal,previous-digest-presence,previous-digest-value,entry-count,ordered-entries,entry-shape-class,row-content-child,executable-build-child,dependency-domain-digest,dependency-artifact-hash,recorded-wall-ns,verdict-tag",
    "excluded_fields=none",
    "consumers=CheckpointReceipt::digest,checkpoint_params,load_chain,checkpoint_staleness_history,staleness_at_checkpointed_with",
    "mutations=digest-domain:crates/fs-roofline/src/checkpoint.rs#staleness_checkpoint_chain_identity_fields_move_independently,identity-version:crates/fs-roofline/src/checkpoint.rs#staleness_checkpoint_chain_identity_versions_fail_closed,checkpoint-body-schema:crates/fs-roofline/src/checkpoint.rs#staleness_checkpoint_chain_identity_fields_move_independently,canonical-json-layout:crates/fs-roofline/src/checkpoint.rs#staleness_checkpoint_chain_identity_fields_move_independently,lowercase-hex-string-framing:crates/fs-roofline/src/checkpoint.rs#checkpoint_chain_string_fields_are_injective_and_canonical,previous-digest-framing:crates/fs-roofline/src/checkpoint.rs#staleness_checkpoint_chain_identity_fields_move_independently,legacy-v1-chain-binding:crates/fs-roofline/src/checkpoint.rs#legacy_v1_chain_is_bound_into_v2_genesis,kernel:crates/fs-roofline/src/checkpoint.rs#staleness_checkpoint_chain_identity_fields_move_independently,kernel-version:crates/fs-roofline/src/checkpoint.rs#staleness_checkpoint_chain_identity_fields_move_independently,checkpoint-ordinal:crates/fs-roofline/src/checkpoint.rs#staleness_checkpoint_chain_identity_fields_move_independently,previous-digest-presence:crates/fs-roofline/src/checkpoint.rs#staleness_checkpoint_chain_identity_fields_move_independently,previous-digest-value:crates/fs-roofline/src/checkpoint.rs#staleness_checkpoint_chain_identity_fields_move_independently,entry-count:crates/fs-roofline/src/checkpoint.rs#staleness_checkpoint_chain_identity_fields_move_independently,ordered-entries:crates/fs-roofline/src/checkpoint.rs#staleness_checkpoint_chain_identity_fields_move_independently,entry-shape-class:crates/fs-roofline/src/checkpoint.rs#staleness_checkpoint_chain_identity_fields_move_independently,row-content-child:crates/fs-roofline/src/checkpoint.rs#staleness_checkpoint_chain_identity_fields_move_independently,executable-build-child:crates/fs-roofline/src/checkpoint.rs#staleness_checkpoint_chain_identity_fields_move_independently,dependency-domain-digest:crates/fs-roofline/src/checkpoint.rs#staleness_checkpoint_chain_identity_fields_move_independently,dependency-artifact-hash:crates/fs-roofline/src/checkpoint.rs#staleness_checkpoint_chain_identity_fields_move_independently,recorded-wall-ns:crates/fs-roofline/src/checkpoint.rs#staleness_checkpoint_chain_identity_fields_move_independently,verdict-tag:crates/fs-roofline/src/checkpoint.rs#staleness_checkpoint_chain_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_staleness_checkpoint_chain_identity_fields",
    "transport_guard=load_chain",
    "version_guard=crates/fs-roofline/src/checkpoint.rs#staleness_checkpoint_chain_identity_versions_fail_closed",
    "coupling_surface=fs-roofline:staleness-checkpoint-chain",
];

/// Per-row verdict sealed into a checkpoint.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CheckpointVerdict {
    /// The row passed the exhaustive verifier at checkpoint time.
    Valid,
    /// The row failed verification — a permanent tombstone.
    Corrupt,
}

impl CheckpointVerdict {
    fn as_str(self) -> &'static str {
        match self {
            CheckpointVerdict::Valid => "valid",
            CheckpointVerdict::Corrupt => "corrupt",
        }
    }

    fn parse(text: &str) -> Option<Self> {
        match text {
            "valid" => Some(CheckpointVerdict::Valid),
            "corrupt" => Some(CheckpointVerdict::Corrupt),
            _ => None,
        }
    }
}

/// One covered row inside a checkpoint body.
#[derive(Debug, Clone, PartialEq, Eq)]
struct CheckpointEntry {
    shape_class: String,
    row_hash: fs_blake3::ContentHash,
    build: fs_blake3::ContentHash,
    dep_digest: fs_blake3::ContentHash,
    dep_artifact: fs_blake3::ContentHash,
    recorded_at_ns: i64,
    verdict: CheckpointVerdict,
}

/// Exact tune-row fields covered by the durable row-content receipt.
struct StalenessRowContentIdentityInput<'a> {
    kernel: &'a [u8],
    shape_class: &'a [u8],
    machine: &'a [u8],
    params: &'a [u8],
    measured: &'a [u8],
}

impl<'a> StalenessRowContentIdentityInput<'a> {
    fn from_row(row: &'a fs_ledger::TuneRow) -> Self {
        Self {
            kernel: row.kernel.as_bytes(),
            shape_class: row.shape_class.as_bytes(),
            machine: row.machine.as_slice(),
            params: row.params.as_bytes(),
            measured: row.measured.as_bytes(),
        }
    }
}

#[allow(dead_code)]
fn classify_staleness_row_content_identity_fields(input: &StalenessRowContentIdentityInput<'_>) {
    let StalenessRowContentIdentityInput {
        kernel,
        shape_class,
        machine,
        params,
        measured,
    } = input;
    let _ = (kernel, shape_class, machine, params, measured);
}

/// Semantic checkpoint body fields covered by the durable chain receipt.
struct StalenessCheckpointChainIdentityInput<'a> {
    kernel: &'a str,
    version: &'a str,
    ordinal: u64,
    previous: Option<fs_blake3::ContentHash>,
    entries: &'a [CheckpointEntry],
}

#[allow(dead_code)]
fn classify_staleness_checkpoint_chain_identity_fields(
    input: &StalenessCheckpointChainIdentityInput<'_>,
) {
    let StalenessCheckpointChainIdentityInput {
        kernel,
        version,
        ordinal,
        previous,
        entries,
    } = input;
    let _ = (kernel, version, ordinal, previous, entries);
    for entry in *entries {
        let CheckpointEntry {
            shape_class,
            row_hash,
            build,
            dep_digest,
            dep_artifact,
            recorded_at_ns,
            verdict,
        } = entry;
        let _ = (
            shape_class,
            row_hash,
            build,
            dep_digest,
            dep_artifact,
            recorded_at_ns,
            verdict,
        );
    }
}

/// A verified checkpoint chain head plus its covered entries.
#[derive(Debug)]
struct VerifiedCheckpoint {
    entries: Vec<CheckpointEntry>,
    /// Digest of the newest checkpoint (the chain tip).
    tip_digest: fs_blake3::ContentHash,
    /// Number of checkpoints in the chain (== next free ordinal).
    len: u64,
}

/// Receipt returned by [`checkpoint_staleness_history`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CheckpointReceipt {
    /// Ordinal of the checkpoint written (0-based, contiguous per chain).
    pub ordinal: u64,
    /// Chained digest sealing this checkpoint.
    pub digest: fs_blake3::ContentHash,
    /// Rows covered.
    pub rows: usize,
    /// Rows sealed as corrupt tombstones (including inherited ones).
    pub corrupt_rows: usize,
}

/// All-zero placeholder hash for tombstone entries (no verified metadata).
fn zero_hash() -> fs_blake3::ContentHash {
    fs_blake3::ContentHash::from_hex(&"0".repeat(64)).expect("zero hash parses")
}

fn checkpoint_kernel(kernel: &str) -> String {
    format!("{CHECKPOINT_KERNEL_PREFIX}{kernel}")
}

/// One sealed entry from an exhaustive-verifier result. `None` (failed or
/// inherited-tombstone) seals a permanent Corrupt tombstone with zeroed
/// metadata.
fn seal_entry(
    row: &fs_ledger::TuneRow,
    validated: Option<ValidatedRooflineRow>,
) -> CheckpointEntry {
    match validated {
        Some(v) => CheckpointEntry {
            shape_class: row.shape_class.clone(),
            row_hash: row_content_hash(row),
            build: v.build_identity,
            dep_digest: v.dependency_receipt_digest,
            dep_artifact: v.dependency_receipt_artifact,
            recorded_at_ns: v.recorded_at_ns,
            verdict: CheckpointVerdict::Valid,
        },
        None => CheckpointEntry {
            shape_class: row.shape_class.clone(),
            row_hash: row_content_hash(row),
            build: zero_hash(),
            dep_digest: zero_hash(),
            dep_artifact: zero_hash(),
            recorded_at_ns: 0,
            verdict: CheckpointVerdict::Corrupt,
        },
    }
}

fn checkpoint_shape_class(version: &str, ordinal: u64) -> String {
    format!("{CHECKPOINT_SHAPE_PREFIX}{version}:{ordinal:08}")
}

fn checkpoint_params(ordinal: u64, digest: fs_blake3::ContentHash) -> String {
    format!(
        "{{\"schema\":\"{CHECKPOINT_SCHEMA}\",\"ordinal\":{ordinal},\"digest\":\"{}\"}}",
        digest.to_hex()
    )
}

fn legacy_checkpoint_shape_class(version: &str, ordinal: u64) -> String {
    format!("{LEGACY_CHECKPOINT_SHAPE_PREFIX}{version}:{ordinal:08}")
}

fn legacy_checkpoint_params(ordinal: u64, digest: fs_blake3::ContentHash) -> String {
    format!(
        "{{\"schema\":\"{LEGACY_CHECKPOINT_SCHEMA}\",\"ordinal\":{ordinal},\"digest\":\"{digest}\"}}"
    )
}

fn push_row_content_field(material: &mut Vec<u8>, field: &[u8]) {
    material.extend_from_slice(&u64::try_from(field.len()).unwrap_or(u64::MAX).to_le_bytes());
    material.extend_from_slice(field);
}

fn staleness_row_content_receipt(
    input: &StalenessRowContentIdentityInput<'_>,
) -> fs_blake3::ContentHash {
    staleness_row_content_receipt_with_domain(ROW_CONTENT_DOMAIN, input)
}

fn staleness_row_content_receipt_with_domain(
    domain: &str,
    input: &StalenessRowContentIdentityInput<'_>,
) -> fs_blake3::ContentHash {
    let mut material = Vec::new();
    for part in [
        input.kernel,
        input.shape_class,
        input.machine,
        input.params,
        input.measured,
    ] {
        push_row_content_field(&mut material, part);
    }
    fs_blake3::hash_domain(domain, &material)
}

/// Content hash binding one tune row's full stored identity.
fn row_content_hash(row: &fs_ledger::TuneRow) -> fs_blake3::ContentHash {
    staleness_row_content_receipt(&StalenessRowContentIdentityInput::from_row(row))
}

fn row_content_is_current(row: &fs_ledger::TuneRow, expected: fs_blake3::ContentHash) -> bool {
    row_content_hash(row) == expected
}

fn lower_hex(bytes: &[u8]) -> String {
    const DIGITS: &[u8; 16] = b"0123456789abcdef";
    let mut encoded = String::with_capacity(bytes.len().saturating_mul(2));
    for byte in bytes {
        encoded.push(char::from(DIGITS[usize::from(byte >> 4)]));
        encoded.push(char::from(DIGITS[usize::from(byte & 0x0f)]));
    }
    encoded
}

fn decode_lower_hex_utf8(encoded: &str) -> Option<String> {
    fn nibble(byte: u8) -> Option<u8> {
        match byte {
            b'0'..=b'9' => Some(byte - b'0'),
            b'a'..=b'f' => Some(byte - b'a' + 10),
            _ => None,
        }
    }

    let bytes = encoded.as_bytes();
    if bytes.len() % 2 != 0 {
        return None;
    }
    let mut decoded = Vec::with_capacity(bytes.len() / 2);
    for pair in bytes.chunks_exact(2) {
        decoded.push((nibble(pair[0])? << 4) | nibble(pair[1])?);
    }
    String::from_utf8(decoded).ok()
}

fn entry_json(entry: &CheckpointEntry) -> String {
    format!(
        "{{\"shape_class_hex\":\"{}\",\"row_hash\":\"{}\",\"build\":\"{}\",\"dep_digest\":\"{}\",\"dep_artifact\":\"{}\",\"recorded_at_ns\":{},\"verdict\":\"{}\"}}",
        lower_hex(entry.shape_class.as_bytes()),
        entry.row_hash.to_hex(),
        entry.build.to_hex(),
        entry.dep_digest.to_hex(),
        entry.dep_artifact.to_hex(),
        entry.recorded_at_ns,
        entry.verdict.as_str(),
    )
}

fn checkpoint_chain_body_json(input: &StalenessCheckpointChainIdentityInput<'_>) -> String {
    checkpoint_chain_body_json_with_schema(CHECKPOINT_SCHEMA, input)
}

fn checkpoint_chain_body_json_with_schema(
    schema: &str,
    input: &StalenessCheckpointChainIdentityInput<'_>,
) -> String {
    let prev_text = input
        .previous
        .map_or_else(|| "null".to_string(), |p| format!("\"{}\"", p.to_hex()));
    let rows = input
        .entries
        .iter()
        .map(entry_json)
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"schema_hex\":\"{}\",\"kernel_hex\":\"{}\",\"version_hex\":\"{}\",\"ordinal\":{},\"prev\":{prev_text},\"rows\":[{rows}]}}",
        lower_hex(schema.as_bytes()),
        lower_hex(input.kernel.as_bytes()),
        lower_hex(input.version.as_bytes()),
        input.ordinal,
    )
}

fn staleness_checkpoint_chain_receipt(
    input: &StalenessCheckpointChainIdentityInput<'_>,
) -> fs_blake3::ContentHash {
    staleness_checkpoint_chain_receipt_with_domain(CHECKPOINT_CHAIN_DOMAIN, input)
}

fn staleness_checkpoint_chain_receipt_with_domain(
    domain: &str,
    input: &StalenessCheckpointChainIdentityInput<'_>,
) -> fs_blake3::ContentHash {
    let body = checkpoint_chain_body_json(input);
    let mut material = Vec::new();
    if let Some(prev) = input.previous {
        material.extend_from_slice(prev.as_bytes());
    }
    material.extend_from_slice(body.as_bytes());
    fs_blake3::hash_domain(domain, &material)
}

fn legacy_entry_json(entry: &CheckpointEntry) -> String {
    format!(
        "{{\"shape_class\":\"{}\",\"row_hash\":\"{}\",\"build\":\"{}\",\"dep_digest\":\"{}\",\"dep_artifact\":\"{}\",\"recorded_at_ns\":{},\"verdict\":\"{}\"}}",
        entry.shape_class,
        entry.row_hash,
        entry.build,
        entry.dep_digest,
        entry.dep_artifact,
        entry.recorded_at_ns,
        entry.verdict.as_str(),
    )
}

fn legacy_body_json(
    kernel: &str,
    version: &str,
    ordinal: u64,
    previous: Option<fs_blake3::ContentHash>,
    entries: &[CheckpointEntry],
) -> String {
    let prev_text = previous.map_or_else(|| "null".to_string(), |p| format!("\"{p}\""));
    let rows = entries
        .iter()
        .map(legacy_entry_json)
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"schema\":\"{LEGACY_CHECKPOINT_SCHEMA}\",\"kernel\":\"{kernel}\",\"version\":\"{version}\",\"ordinal\":{ordinal},\"prev\":{prev_text},\"rows\":[{rows}]}}",
    )
}

fn legacy_chain_digest(
    previous: Option<fs_blake3::ContentHash>,
    body: &str,
) -> fs_blake3::ContentHash {
    let mut material = Vec::new();
    if let Some(previous) = previous {
        material.extend_from_slice(previous.as_bytes());
    }
    material.extend_from_slice(body.as_bytes());
    fs_blake3::hash_domain(LEGACY_CHECKPOINT_CHAIN_DOMAIN, &material)
}

/// Parsed checkpoint body fields.
struct ParsedBody {
    kernel: String,
    version: String,
    ordinal: u64,
    prev: Option<fs_blake3::ContentHash>,
    entries: Vec<CheckpointEntry>,
}

/// Strict parse of one checkpoint body; byte-exact round trip enforced.
fn parse_body(text: &str) -> Option<ParsedBody> {
    let rest = text.strip_prefix(&format!(
        "{{\"schema_hex\":\"{}\",\"kernel_hex\":\"",
        lower_hex(CHECKPOINT_SCHEMA.as_bytes())
    ))?;
    let (kernel_hex, rest) = rest.split_once("\",\"version_hex\":\"")?;
    let (version_hex, rest) = rest.split_once("\",\"ordinal\":")?;
    let kernel = decode_lower_hex_utf8(kernel_hex)?;
    let version = decode_lower_hex_utf8(version_hex)?;
    let (ordinal_text, rest) = rest.split_once(",\"prev\":")?;
    let ordinal: u64 = ordinal_text.parse().ok()?;
    let (prev_text, rest) = rest.split_once(",\"rows\":[")?;
    let prev = if prev_text == "null" {
        None
    } else {
        Some(fs_blake3::ContentHash::from_hex(
            prev_text.strip_prefix('"')?.strip_suffix('"')?,
        )?)
    };
    let rows_text = rest.strip_suffix("]}")?;
    let mut entries = Vec::new();
    if !rows_text.is_empty() {
        for raw in rows_text.split("},") {
            let raw = if raw.ends_with('}') {
                raw.to_string()
            } else {
                format!("{raw}}}")
            };
            let inner = raw.strip_prefix("{\"shape_class_hex\":\"")?;
            let (shape_class_hex, inner) = inner.split_once("\",\"row_hash\":\"")?;
            let (row_hash, inner) = inner.split_once("\",\"build\":\"")?;
            let (build, inner) = inner.split_once("\",\"dep_digest\":\"")?;
            let (dep_digest, inner) = inner.split_once("\",\"dep_artifact\":\"")?;
            let (dep_artifact, inner) = inner.split_once("\",\"recorded_at_ns\":")?;
            let (recorded_text, inner) = inner.split_once(",\"verdict\":\"")?;
            let verdict_text = inner.strip_suffix("\"}")?;
            entries.push(CheckpointEntry {
                shape_class: decode_lower_hex_utf8(shape_class_hex)?,
                row_hash: fs_blake3::ContentHash::from_hex(row_hash)?,
                build: fs_blake3::ContentHash::from_hex(build)?,
                dep_digest: fs_blake3::ContentHash::from_hex(dep_digest)?,
                dep_artifact: fs_blake3::ContentHash::from_hex(dep_artifact)?,
                recorded_at_ns: recorded_text.parse().ok()?,
                verdict: CheckpointVerdict::parse(verdict_text)?,
            });
        }
    }
    // Byte-exact round trip: non-canonical spellings are refused.
    let identity_input = StalenessCheckpointChainIdentityInput {
        kernel: &kernel,
        version: &version,
        ordinal,
        previous: prev,
        entries: &entries,
    };
    (checkpoint_chain_body_json(&identity_input) == text).then_some(ParsedBody {
        kernel,
        version,
        ordinal,
        prev,
        entries,
    })
}

/// Frozen decoder for the v1 checkpoint body. It intentionally preserves the
/// original raw-string layout; v2 never rewrites legacy bytes under new rules.
fn parse_legacy_body(text: &str) -> Option<ParsedBody> {
    let rest = text.strip_prefix(&format!(
        "{{\"schema\":\"{LEGACY_CHECKPOINT_SCHEMA}\",\"kernel\":\""
    ))?;
    let (kernel, rest) = rest.split_once("\",\"version\":\"")?;
    let (version, rest) = rest.split_once("\",\"ordinal\":")?;
    let (ordinal_text, rest) = rest.split_once(",\"prev\":")?;
    let ordinal: u64 = ordinal_text.parse().ok()?;
    let (prev_text, rest) = rest.split_once(",\"rows\":[")?;
    let prev = if prev_text == "null" {
        None
    } else {
        Some(fs_blake3::ContentHash::from_hex(
            prev_text.strip_prefix('"')?.strip_suffix('"')?,
        )?)
    };
    let rows_text = rest.strip_suffix("]}")?;
    let mut entries = Vec::new();
    if !rows_text.is_empty() {
        for raw in rows_text.split("},") {
            let raw = if raw.ends_with('}') {
                raw.to_string()
            } else {
                format!("{raw}}}")
            };
            let inner = raw.strip_prefix("{\"shape_class\":\"")?;
            let (shape_class, inner) = inner.split_once("\",\"row_hash\":\"")?;
            let (row_hash, inner) = inner.split_once("\",\"build\":\"")?;
            let (build, inner) = inner.split_once("\",\"dep_digest\":\"")?;
            let (dep_digest, inner) = inner.split_once("\",\"dep_artifact\":\"")?;
            let (dep_artifact, inner) = inner.split_once("\",\"recorded_at_ns\":")?;
            let (recorded_text, inner) = inner.split_once(",\"verdict\":\"")?;
            let verdict_text = inner.strip_suffix("\"}")?;
            entries.push(CheckpointEntry {
                shape_class: shape_class.to_string(),
                row_hash: fs_blake3::ContentHash::from_hex(row_hash)?,
                build: fs_blake3::ContentHash::from_hex(build)?,
                dep_digest: fs_blake3::ContentHash::from_hex(dep_digest)?,
                dep_artifact: fs_blake3::ContentHash::from_hex(dep_artifact)?,
                recorded_at_ns: recorded_text.parse().ok()?,
                verdict: CheckpointVerdict::parse(verdict_text)?,
            });
        }
    }
    (legacy_body_json(kernel, version, ordinal, prev, &entries) == text).then_some(ParsedBody {
        kernel: kernel.to_string(),
        version: version.to_string(),
        ordinal,
        prev,
        entries,
    })
}

/// Chain-load outcome: absent, present-but-unverifiable, or verified.
enum ChainState {
    /// No checkpoint rows exist for this (kernel, version, machine).
    Empty,
    /// Rows exist but the chain fails verification (parse failure, ordinal
    /// gap, digest mismatch, params/shape inconsistency).
    Broken,
    /// The full chain verified; carries the newest checkpoint's entries.
    Verified(VerifiedCheckpoint),
}

fn checkpoint_ordinal(shape_class: &str, prefix: &str, version: &str) -> Result<Option<u64>, ()> {
    let Some(suffix) = shape_class.strip_prefix(prefix) else {
        return Ok(None);
    };
    let Some((stored_version, ordinal)) = suffix.rsplit_once(':') else {
        return if suffix == version { Err(()) } else { Ok(None) };
    };
    if stored_version != version {
        return Ok(None);
    }
    if ordinal.len() != 8 || !ordinal.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err(());
    }
    ordinal.parse().map(Some).map_err(|_| ())
}

fn entries_are_canonical(entries: &[CheckpointEntry]) -> bool {
    !entries.is_empty()
        && entries
            .windows(2)
            .all(|pair| pair[0].shape_class < pair[1].shape_class)
}

fn corrupt_entry_is_canonical(entry: &CheckpointEntry) -> bool {
    entry.verdict == CheckpointVerdict::Corrupt
        && entry.build == zero_hash()
        && entry.dep_digest == zero_hash()
        && entry.dep_artifact == zero_hash()
        && entry.recorded_at_ns == 0
}

fn canonical_tombstone(shape_class: String, row_hash: fs_blake3::ContentHash) -> CheckpointEntry {
    CheckpointEntry {
        shape_class,
        row_hash,
        build: zero_hash(),
        dep_digest: zero_hash(),
        dep_artifact: zero_hash(),
        recorded_at_ns: 0,
        verdict: CheckpointVerdict::Corrupt,
    }
}

fn entry_transition_is_monotonic(previous: &CheckpointEntry, current: &CheckpointEntry) -> bool {
    if previous.shape_class != current.shape_class {
        return false;
    }
    match (previous.verdict, current.verdict) {
        (CheckpointVerdict::Valid, CheckpointVerdict::Valid)
        | (CheckpointVerdict::Corrupt, CheckpointVerdict::Corrupt) => previous == current,
        (CheckpointVerdict::Valid, CheckpointVerdict::Corrupt) => {
            current.row_hash == previous.row_hash && corrupt_entry_is_canonical(current)
        }
        (CheckpointVerdict::Corrupt, CheckpointVerdict::Valid) => false,
    }
}

fn entries_extend_monotonically(previous: &[CheckpointEntry], current: &[CheckpointEntry]) -> bool {
    previous.iter().all(|old| {
        current
            .binary_search_by(|candidate| candidate.shape_class.cmp(&old.shape_class))
            .ok()
            .is_some_and(|index| entry_transition_is_monotonic(old, &current[index]))
    })
}

/// Fold one historically valid v1 snapshot into a monotonic migration view.
///
/// The v1 writer authenticated every raw snapshot but did not authenticate
/// transitions between snapshots. In particular, it recomputed a tombstone's
/// row hash from whatever bytes were present at each reseal, and a row absent
/// from one snapshot could later reappear. Rejecting those histories would
/// strand genuine retained v1 evidence. Instead, v2 binds the exact v1 tip and
/// carries this conservative union forward: every shape ever seen remains
/// covered, only byte-identical valid metadata remains valid, and the first
/// corrupt, rewritten, or omitted verdict becomes a permanent canonical
/// tombstone whose row hash never moves again.
fn normalize_legacy_entries(
    normalized: &mut std::collections::BTreeMap<String, CheckpointEntry>,
    snapshot: &[CheckpointEntry],
) -> bool {
    let mut observed = std::collections::BTreeSet::new();
    for current in snapshot {
        if current.verdict == CheckpointVerdict::Corrupt && !corrupt_entry_is_canonical(current) {
            return false;
        }
        observed.insert(current.shape_class.clone());
        match normalized.get(&current.shape_class).cloned() {
            Some(previous) if previous.verdict == CheckpointVerdict::Corrupt => {
                // A v1 restoration could change the stored tombstone hash or
                // even reappear as valid after an intervening absent snapshot.
                // Neither transition is allowed to erase the first incident.
            }
            Some(previous)
                if current.verdict == CheckpointVerdict::Corrupt || previous != *current =>
            {
                normalized.insert(
                    current.shape_class.clone(),
                    canonical_tombstone(current.shape_class.clone(), previous.row_hash),
                );
            }
            _ => {
                normalized.insert(current.shape_class.clone(), current.clone());
            }
        }
    }
    for (shape_class, previous) in normalized.iter_mut() {
        if previous.verdict == CheckpointVerdict::Valid && !observed.contains(shape_class) {
            *previous = canonical_tombstone(shape_class.clone(), previous.row_hash);
        }
    }
    true
}

fn verify_legacy_chain(
    rows: &[fs_ledger::TuneRow],
    kernel: &str,
    version: &str,
    machine_key: [u8; 40],
) -> ChainState {
    let mut chain = Vec::new();
    for row in rows.iter().filter(|row| row.machine == machine_key) {
        match checkpoint_ordinal(&row.shape_class, LEGACY_CHECKPOINT_SHAPE_PREFIX, version) {
            Ok(Some(ordinal)) => chain.push((ordinal, row)),
            Ok(None) => {}
            Err(()) => return ChainState::Broken,
        }
    }
    if chain.is_empty() {
        return ChainState::Empty;
    }
    chain.sort_by_key(|(ordinal, _)| *ordinal);
    let mut previous_digest = None;
    let mut normalized_entries = std::collections::BTreeMap::new();
    for (index, (shape_ordinal, row)) in chain.iter().enumerate() {
        let Some(body) = parse_legacy_body(&row.measured) else {
            return ChainState::Broken;
        };
        let expected_digest = legacy_chain_digest(previous_digest, &row.measured);
        if *shape_ordinal != index as u64
            || body.kernel != kernel
            || body.version != version
            || body.ordinal != index as u64
            || body.prev != previous_digest
            || row.params != legacy_checkpoint_params(body.ordinal, expected_digest)
            || row.shape_class != legacy_checkpoint_shape_class(version, body.ordinal)
            || !entries_are_canonical(&body.entries)
            || !normalize_legacy_entries(&mut normalized_entries, &body.entries)
        {
            return ChainState::Broken;
        }
        previous_digest = Some(expected_digest);
    }
    match previous_digest {
        Some(tip_digest) if !normalized_entries.is_empty() => {
            ChainState::Verified(VerifiedCheckpoint {
                entries: normalized_entries.into_values().collect(),
                tip_digest,
                len: chain.len() as u64,
            })
        }
        _ => ChainState::Broken,
    }
}

fn verify_current_chain(
    rows: &[fs_ledger::TuneRow],
    kernel: &str,
    version: &str,
    machine_key: [u8; 40],
    legacy_tip: Option<fs_blake3::ContentHash>,
    legacy_entries: Option<&[CheckpointEntry]>,
) -> ChainState {
    let mut chain = Vec::new();
    for row in rows.iter().filter(|row| row.machine == machine_key) {
        match checkpoint_ordinal(&row.shape_class, CHECKPOINT_SHAPE_PREFIX, version) {
            Ok(Some(ordinal)) => chain.push((ordinal, row)),
            Ok(None) => {}
            Err(()) => return ChainState::Broken,
        }
    }
    if chain.is_empty() {
        return ChainState::Empty;
    }
    chain.sort_by_key(|(ordinal, _)| *ordinal);
    let mut previous_digest = legacy_tip;
    let mut previous_entries = legacy_entries.map(|entries| entries.to_vec());
    for (index, (shape_ordinal, row)) in chain.iter().enumerate() {
        let Some(body) = parse_body(&row.measured) else {
            return ChainState::Broken;
        };
        let expected_digest =
            staleness_checkpoint_chain_receipt(&StalenessCheckpointChainIdentityInput {
                kernel: &body.kernel,
                version: &body.version,
                ordinal: body.ordinal,
                previous: body.prev,
                entries: &body.entries,
            });
        if *shape_ordinal != index as u64
            || body.kernel != kernel
            || body.version != version
            || body.ordinal != index as u64
            || body.prev != previous_digest
            || row.params != checkpoint_params(body.ordinal, expected_digest)
            || row.shape_class != checkpoint_shape_class(version, body.ordinal)
            || !entries_are_canonical(&body.entries)
            || previous_entries
                .as_deref()
                .is_some_and(|prior| !entries_extend_monotonically(prior, &body.entries))
        {
            return ChainState::Broken;
        }
        previous_digest = Some(expected_digest);
        previous_entries = Some(body.entries);
    }
    match (previous_entries, previous_digest) {
        (Some(entries), Some(tip_digest)) => ChainState::Verified(VerifiedCheckpoint {
            entries,
            tip_digest,
            len: chain.len() as u64,
        }),
        _ => ChainState::Broken,
    }
}

/// Load and verify the checkpoint chain for `(kernel, version, machine)`.
/// Callers FAIL CLOSED on both `Empty` and `Broken`.
fn load_chain(
    ledger: &Ledger,
    kernel: &str,
    version: &str,
    machine_key: [u8; 40],
) -> Result<ChainState, LedgerError> {
    let rows = ledger.tune_rows(&checkpoint_kernel(kernel))?;
    let legacy = verify_legacy_chain(&rows, kernel, version, machine_key);
    let (legacy_tip, legacy_entries) = match legacy {
        ChainState::Empty => (None, None),
        ChainState::Broken => return Ok(ChainState::Broken),
        ChainState::Verified(verified) => (Some(verified.tip_digest), Some(verified.entries)),
    };
    match verify_current_chain(
        &rows,
        kernel,
        version,
        machine_key,
        legacy_tip,
        legacy_entries.as_deref(),
    ) {
        ChainState::Verified(verified) => Ok(ChainState::Verified(verified)),
        ChainState::Broken => Ok(ChainState::Broken),
        ChainState::Empty => match (legacy_tip, legacy_entries) {
            (Some(tip_digest), Some(entries)) => {
                Ok(ChainState::Verified(VerifiedCheckpoint {
                    entries,
                    tip_digest,
                    // v2 starts from ordinal zero even when its genesis is
                    // chained to a verified v1 tip.
                    len: 0,
                }))
            }
            _ => Ok(ChainState::Empty),
        },
    }
}

/// Run the exhaustive verifier once over the current matching rows and seal
/// the results as the next checkpoint in the chain. Prior tombstones are
/// PRESERVED: a row sealed corrupt in any earlier checkpoint remains corrupt
/// here even if its stored bytes now re-verify (never un-corrupt).
///
/// # Errors
/// Ledger errors propagate. A lattice-prefix verdict (no rows for this
/// machine/baseline) is an [`LedgerError::Invalid`] refusal — there is
/// nothing to seal.
pub fn checkpoint_staleness_history(
    ledger: &Ledger,
    kernel: &str,
    version: &str,
    current_fingerprint: u64,
    current_baseline: fs_blake3::ContentHash,
) -> Result<CheckpointReceipt, LedgerError> {
    checkpoint_history_with_dependency(
        ledger,
        kernel,
        version,
        current_fingerprint,
        current_baseline,
        DependencyReceiptBinding::current().ok(),
    )
}

pub(crate) fn checkpoint_history_with_dependency(
    ledger: &Ledger,
    kernel: &str,
    version: &str,
    current_fingerprint: u64,
    current_baseline: fs_blake3::ContentHash,
    expected_dependency: Option<DependencyReceiptBinding>,
) -> Result<CheckpointReceipt, LedgerError> {
    if ledger.in_transaction() {
        return Err(LedgerError::Invalid {
            field: "staleness-checkpoint-transaction".to_string(),
            problem: "checkpoint sealing owns its rollback boundary; commit or roll back the \
                      caller transaction before sealing"
                .to_string(),
        });
    }
    ledger.begin()?;
    let result = checkpoint_history_with_dependency_inner(
        ledger,
        kernel,
        version,
        current_fingerprint,
        current_baseline,
        expected_dependency,
    );
    match result {
        Ok(receipt) => match ledger.commit() {
            Ok(()) => Ok(receipt),
            Err(error) => Err(rollback_checkpoint_transaction(ledger, error)),
        },
        Err(error) => Err(rollback_checkpoint_transaction(ledger, error)),
    }
}

fn rollback_checkpoint_transaction(ledger: &Ledger, primary: LedgerError) -> LedgerError {
    match ledger.rollback() {
        Ok(()) => primary,
        Err(rollback) => LedgerError::Invalid {
            field: "staleness-checkpoint-transaction".to_string(),
            problem: format!(
                "checkpoint operation failed ({primary}); rollback also failed ({rollback}); \
                 transaction state is uncertain"
            ),
        },
    }
}

fn checkpoint_history_with_dependency_inner(
    ledger: &Ledger,
    kernel: &str,
    version: &str,
    current_fingerprint: u64,
    current_baseline: fs_blake3::ContentHash,
    expected_dependency: Option<DependencyReceiptBinding>,
) -> Result<CheckpointReceipt, LedgerError> {
    let machine_key = roofline_machine_key(current_fingerprint, current_baseline);
    // Load history before the production-row selection. That order prevents
    // total rollback from hiding a verified chain behind NeverMeasured.
    let (prior_entries, next_ordinal, prev_digest) =
        match load_chain(ledger, kernel, version, machine_key)? {
            ChainState::Verified(verified) => {
                (verified.entries, verified.len, Some(verified.tip_digest))
            }
            ChainState::Empty => (Vec::new(), 0, None),
            ChainState::Broken => {
                return Err(LedgerError::Invalid {
                    field: "staleness-checkpoint".to_string(),
                    problem: "existing checkpoint history fails verification; refusing to \
                              extend it (broken v1/v2 continuity is permanent evidence, \
                              never overwritten)"
                        .to_string(),
                });
            }
        };
    let matching = match select_matching_rows(
        ledger,
        kernel,
        version,
        current_fingerprint,
        Some(current_baseline),
    )? {
        RowSelection::Rows(rows) => rows,
        RowSelection::Verdict(_verdict)
            if !prior_entries.is_empty()
                && prior_entries
                    .iter()
                    .all(|entry| entry.verdict == CheckpointVerdict::Corrupt) =>
        {
            Vec::new()
        }
        RowSelection::Verdict(verdict) => {
            return Err(LedgerError::Invalid {
                field: "staleness-checkpoint".to_string(),
                problem: format!(
                    "cannot seal after row selection classified {verdict:?}; verified valid \
                     history must not disappear"
                ),
            });
        }
    };

    let mut prior_by_shape: std::collections::BTreeMap<String, CheckpointEntry> = prior_entries
        .into_iter()
        .map(|entry| (entry.shape_class.clone(), entry))
        .collect();
    let mut entries = Vec::with_capacity(matching.len().max(prior_by_shape.len()));
    let mut corrupt_rows = 0usize;
    for row in &matching {
        if let Some(previous) = prior_by_shape.remove(&row.shape_class) {
            if previous.verdict == CheckpointVerdict::Corrupt {
                corrupt_rows += 1;
                entries.push(previous);
                continue;
            }
            if !row_content_is_current(row, previous.row_hash) {
                return Err(LedgerError::Invalid {
                    field: "staleness-checkpoint".to_string(),
                    problem: format!(
                        "previously valid covered row {:?} changed; refusing to launder \
                         mutation into a newer checkpoint",
                        row.shape_class
                    ),
                });
            }
            let validated = validate_roofline_row(
                ledger,
                row,
                kernel,
                version,
                current_fingerprint,
                current_baseline,
                expected_dependency,
            )?;
            let entry = seal_entry(row, validated);
            if entry.verdict == CheckpointVerdict::Corrupt {
                corrupt_rows += 1;
            } else if entry != previous {
                return Err(LedgerError::Invalid {
                    field: "staleness-checkpoint".to_string(),
                    problem: format!(
                        "previously valid covered row {:?} revalidated to different metadata",
                        row.shape_class
                    ),
                });
            }
            entries.push(entry);
            continue;
        }
        let validated = validate_roofline_row(
            ledger,
            row,
            kernel,
            version,
            current_fingerprint,
            current_baseline,
            expected_dependency,
        )?;
        if validated.is_none() {
            corrupt_rows += 1;
        }
        entries.push(seal_entry(row, validated));
    }
    for (_, previous) in prior_by_shape {
        if previous.verdict == CheckpointVerdict::Valid {
            return Err(LedgerError::Invalid {
                field: "staleness-checkpoint".to_string(),
                problem: format!(
                    "previously valid covered row {:?} is missing; refusing to reseal \
                     truncated history",
                    previous.shape_class
                ),
            });
        }
        corrupt_rows += 1;
        entries.push(previous);
    }
    entries.sort_by(|a, b| a.shape_class.cmp(&b.shape_class));

    let checkpoint_input = StalenessCheckpointChainIdentityInput {
        kernel,
        version,
        ordinal: next_ordinal,
        previous: prev_digest,
        entries: &entries,
    };
    let body = checkpoint_chain_body_json(&checkpoint_input);
    let digest = staleness_checkpoint_chain_receipt(&checkpoint_input);
    let params = checkpoint_params(next_ordinal, digest);
    // Append-only: a colliding ordinal must never overwrite sealed history.
    ledger.tune_put_if_absent(
        &checkpoint_kernel(kernel),
        &checkpoint_shape_class(version, next_ordinal),
        &machine_key,
        &params,
        &body,
    )?;
    let stored = ledger
        .tune_get(
            &checkpoint_kernel(kernel),
            &checkpoint_shape_class(version, next_ordinal),
            &machine_key,
        )?
        .ok_or_else(|| LedgerError::Invalid {
            field: "staleness-checkpoint".to_string(),
            problem: "checkpoint row missing immediately after insert".to_string(),
        })?;
    if stored.measured != body || stored.params != params {
        return Err(LedgerError::Invalid {
            field: "staleness-checkpoint".to_string(),
            problem: format!(
                "checkpoint ordinal {next_ordinal} already sealed with different content; \
                 refusing to overwrite chain history"
            ),
        });
    }
    match load_chain(ledger, kernel, version, machine_key)? {
        ChainState::Verified(verified)
            if verified.tip_digest == digest
                && verified.len == next_ordinal + 1
                && verified.entries == entries =>
        {
            // Strict in-transaction re-read closes the v1-tip/v2-genesis race
            // and checks the exact stored history before commit.
        }
        _ => {
            return Err(LedgerError::Invalid {
                field: "staleness-checkpoint".to_string(),
                problem: "checkpoint history failed strict verification immediately after insert"
                    .to_string(),
            });
        }
    }
    Ok(CheckpointReceipt {
        ordinal: next_ordinal,
        digest,
        rows: entries.len(),
        corrupt_rows,
    })
}

/// Checkpoint-accelerated staleness: identical verdict lattice to
/// [`crate::staleness_at`], with covered rows verified against the sealed
/// chain instead of re-walking their manifests. Falls back to the exhaustive
/// path whenever no verified chain exists (fail closed).
///
/// # Errors
/// Ledger and executable-identity errors propagate.
pub fn staleness_at_checkpointed(
    ledger: &Ledger,
    kernel: &str,
    version: &str,
    current_fingerprint: u64,
    current_baseline: Option<fs_blake3::ContentHash>,
    observed_wall_ns: i64,
) -> Result<Staleness, LedgerError> {
    let current_build = executable_build_identity()?;
    staleness_at_checkpointed_with(
        ledger,
        kernel,
        version,
        current_fingerprint,
        current_baseline,
        observed_wall_ns,
        current_build,
        DependencyReceiptBinding::current().ok(),
    )
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn staleness_at_checkpointed_with(
    ledger: &Ledger,
    kernel: &str,
    version: &str,
    current_fingerprint: u64,
    current_baseline: Option<fs_blake3::ContentHash>,
    observed_wall_ns: i64,
    current_build: fs_blake3::ContentHash,
    expected_dependency: Option<DependencyReceiptBinding>,
) -> Result<Staleness, LedgerError> {
    let chain_state = current_baseline
        .map(|baseline| {
            load_chain(
                ledger,
                kernel,
                version,
                roofline_machine_key(current_fingerprint, baseline),
            )
        })
        .transpose()?;
    let matching = match select_matching_rows(
        ledger,
        kernel,
        version,
        current_fingerprint,
        current_baseline,
    )? {
        RowSelection::Verdict(verdict) => {
            if matches!(
                chain_state.as_ref(),
                Some(ChainState::Verified(verified)) if !verified.entries.is_empty()
            ) {
                return Ok(Staleness::CorruptEvidence);
            }
            return Ok(verdict);
        }
        RowSelection::Rows(rows) => rows,
    };
    let baseline = current_baseline.expect("baseline present when rows match");

    let verified = match chain_state.expect("baseline present when rows match") {
        ChainState::Verified(v) => v,
        // No chain, or a chain that fails verification: FAIL CLOSED to the
        // exhaustive per-row path over the already-selected rows.
        ChainState::Empty | ChainState::Broken => {
            return exhaustive_scan(
                ledger,
                &matching,
                kernel,
                version,
                current_fingerprint,
                baseline,
                observed_wall_ns,
                current_build,
                expected_dependency,
            );
        }
    };

    let mut build_scan = BuildRowScan::default();
    let mut covered: std::collections::BTreeMap<&str, &CheckpointEntry> = verified
        .entries
        .iter()
        .map(|e| (e.shape_class.as_str(), e))
        .collect();
    for row in &matching {
        if let Some(entry) = covered.remove(row.shape_class.as_str()) {
            if entry.verdict == CheckpointVerdict::Corrupt {
                // Tombstone: permanently corrupt.
                return Ok(Staleness::CorruptEvidence);
            }
            if !row_content_is_current(row, entry.row_hash) {
                // Tampered since sealing.
                return Ok(Staleness::CorruptEvidence);
            }
            let dependency_matches_current = expected_dependency.is_some_and(|expected| {
                entry.dep_artifact == expected.artifact_hash
                    && entry.dep_digest == expected.domain_digest
            });
            let replayed = ValidatedRooflineRow {
                build_identity: entry.build,
                recorded_at_ns: entry.recorded_at_ns,
                dependency_matches_current,
                dependency_receipt_digest: entry.dep_digest,
                dependency_receipt_artifact: entry.dep_artifact,
            };
            if !build_scan.observe(&replayed, current_build) {
                return Ok(Staleness::CorruptEvidence);
            }
        } else {
            // Delta row (newer than the checkpoint): full validation.
            let Some(validated) = validate_roofline_row(
                ledger,
                row,
                kernel,
                version,
                current_fingerprint,
                baseline,
                expected_dependency,
            )?
            else {
                return Ok(Staleness::CorruptEvidence);
            };
            if !build_scan.observe(&validated, current_build) {
                return Ok(Staleness::CorruptEvidence);
            }
        }
    }
    if !covered.is_empty() {
        // A sealed row vanished from the ledger: rollback/removal.
        return Ok(Staleness::CorruptEvidence);
    }
    Ok(classify_scanned_rows(build_scan, observed_wall_ns))
}

/// The exhaustive per-row scan over pre-selected rows (shared shape with
/// `staleness_at_with_build_and_dependency`'s tail, reused by the fallback).
#[allow(clippy::too_many_arguments)]
fn exhaustive_scan(
    ledger: &Ledger,
    matching: &[fs_ledger::TuneRow],
    kernel: &str,
    version: &str,
    current_fingerprint: u64,
    baseline: fs_blake3::ContentHash,
    observed_wall_ns: i64,
    current_build: fs_blake3::ContentHash,
    expected_dependency: Option<DependencyReceiptBinding>,
) -> Result<Staleness, LedgerError> {
    let mut build_scan = BuildRowScan::default();
    for row in matching {
        let Some(validated) = validate_roofline_row(
            ledger,
            row,
            kernel,
            version,
            current_fingerprint,
            baseline,
            expected_dependency,
        )?
        else {
            return Ok(Staleness::CorruptEvidence);
        };
        if !build_scan.observe(&validated, current_build) {
            return Ok(Staleness::CorruptEvidence);
        }
    }
    Ok(classify_scanned_rows(build_scan, observed_wall_ns))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernels::default_registry;
    use crate::production::{ProductionProbe, ProductionRunConfig};
    use crate::{
        AttestedAxisBaselinePolicy, BaselineAxes, BaselineCandidate, BaselineIdentity, MachineAxes,
        STALENESS_MAX_AGE_NS, promote_baseline,
    };

    const FINGERPRINT: u64 = 0xC4EC;

    fn row_identity_input() -> StalenessRowContentIdentityInput<'static> {
        StalenessRowContentIdentityInput {
            kernel: b"kernel-a",
            shape_class: b"shape-a",
            machine: b"machine-a",
            params: b"params-a",
            measured: b"measured-a",
        }
    }

    fn checkpoint_identity_entry(tag: u8) -> CheckpointEntry {
        CheckpointEntry {
            shape_class: format!("shape-{tag}"),
            row_hash: fs_blake3::ContentHash([tag; 32]),
            build: fs_blake3::ContentHash([tag.wrapping_add(1); 32]),
            dep_digest: fs_blake3::ContentHash([tag.wrapping_add(2); 32]),
            dep_artifact: fs_blake3::ContentHash([tag.wrapping_add(3); 32]),
            recorded_at_ns: i64::from(tag),
            verdict: CheckpointVerdict::Valid,
        }
    }

    fn checkpoint_identity_receipt(
        kernel: &str,
        version: &str,
        ordinal: u64,
        previous: Option<fs_blake3::ContentHash>,
        entries: &[CheckpointEntry],
    ) -> fs_blake3::ContentHash {
        staleness_checkpoint_chain_receipt(&StalenessCheckpointChainIdentityInput {
            kernel,
            version,
            ordinal,
            previous,
            entries,
        })
    }

    #[test]
    fn staleness_row_content_identity_fields_move_independently() {
        let input = row_identity_input();
        let current = staleness_row_content_receipt(&input);
        let changed_domain = staleness_row_content_receipt_with_domain(
            "org.frankensim.fs-roofline.staleness-row-content-shadow.v1",
            &input,
        );
        assert_ne!(current, changed_domain, "the digest domain is semantic");

        for changed in [
            StalenessRowContentIdentityInput {
                kernel: b"kernel-b",
                ..row_identity_input()
            },
            StalenessRowContentIdentityInput {
                shape_class: b"shape-b",
                ..row_identity_input()
            },
            StalenessRowContentIdentityInput {
                machine: b"machine-b",
                ..row_identity_input()
            },
            StalenessRowContentIdentityInput {
                params: b"params-b",
                ..row_identity_input()
            },
            StalenessRowContentIdentityInput {
                measured: b"measured-b",
                ..row_identity_input()
            },
        ] {
            assert_ne!(
                current,
                staleness_row_content_receipt(&changed),
                "every stored tune-row field is semantic"
            );
        }

        let mut shadow_material = Vec::new();
        for part in [
            input.kernel,
            input.shape_class,
            input.machine,
            input.params,
            input.measured,
        ] {
            shadow_material.extend_from_slice(
                &u32::try_from(part.len())
                    .expect("small identity fixture")
                    .to_le_bytes(),
            );
            shadow_material.extend_from_slice(part);
        }
        assert_ne!(
            current,
            fs_blake3::hash_domain(ROW_CONTENT_DOMAIN, &shadow_material),
            "the u64 length-prefix framing is semantic"
        );
    }

    #[test]
    fn staleness_row_content_identity_versions_fail_closed() {
        assert_eq!(STALENESS_ROW_CONTENT_IDENTITY_VERSION, 1);
        assert!(ROW_CONTENT_DOMAIN.ends_with(".v1"));
        let input = row_identity_input();
        assert_ne!(
            staleness_row_content_receipt(&input),
            staleness_row_content_receipt_with_domain(
                "org.frankensim.fs-roofline.staleness-row-content.v2",
                &input,
            )
        );
    }

    #[test]
    fn staleness_checkpoint_chain_identity_fields_move_independently() {
        let previous = Some(fs_blake3::ContentHash([9; 32]));
        let entries = vec![checkpoint_identity_entry(1), checkpoint_identity_entry(2)];
        let current = checkpoint_identity_receipt("kernel-a", "v1", 7, previous, &entries);
        let input = StalenessCheckpointChainIdentityInput {
            kernel: "kernel-a",
            version: "v1",
            ordinal: 7,
            previous,
            entries: &entries,
        };
        assert_ne!(
            current,
            staleness_checkpoint_chain_receipt_with_domain(
                "org.frankensim.fs-roofline.staleness-checkpoint-chain-shadow.v2",
                &input,
            ),
            "the digest domain is semantic"
        );
        assert_ne!(
            current,
            checkpoint_identity_receipt("kernel-b", "v1", 7, previous, &entries),
            "kernel is semantic"
        );
        assert_ne!(
            current,
            checkpoint_identity_receipt("kernel-a", "v2", 7, previous, &entries),
            "kernel version is semantic"
        );
        assert_ne!(
            current,
            checkpoint_identity_receipt("kernel-a", "v1", 8, previous, &entries),
            "ordinal is semantic"
        );
        assert_ne!(
            current,
            checkpoint_identity_receipt("kernel-a", "v1", 7, None, &entries),
            "previous-digest presence is semantic"
        );
        assert_ne!(
            current,
            checkpoint_identity_receipt(
                "kernel-a",
                "v1",
                7,
                Some(fs_blake3::ContentHash([10; 32])),
                &entries,
            ),
            "previous-digest value is semantic"
        );
        assert_ne!(
            current,
            checkpoint_identity_receipt("kernel-a", "v1", 7, previous, &entries[..1]),
            "entry count is semantic"
        );
        let mut reversed = entries.clone();
        reversed.reverse();
        assert_ne!(
            current,
            checkpoint_identity_receipt("kernel-a", "v1", 7, previous, &reversed),
            "entry order is semantic"
        );

        for mutate in 0..7 {
            let mut changed = entries.clone();
            match mutate {
                0 => changed[0].shape_class.push_str("-changed"),
                1 => changed[0].row_hash = fs_blake3::ContentHash([11; 32]),
                2 => changed[0].build = fs_blake3::ContentHash([12; 32]),
                3 => changed[0].dep_digest = fs_blake3::ContentHash([13; 32]),
                4 => changed[0].dep_artifact = fs_blake3::ContentHash([14; 32]),
                5 => changed[0].recorded_at_ns += 1,
                6 => changed[0].verdict = CheckpointVerdict::Corrupt,
                _ => unreachable!(),
            }
            assert_ne!(
                current,
                checkpoint_identity_receipt("kernel-a", "v1", 7, previous, &changed),
                "every checkpoint entry field is semantic (mutation {mutate})"
            );
        }

        let changed_schema_body = checkpoint_chain_body_json_with_schema(
            "fs-roofline-staleness-checkpoint-shadow-v2",
            &input,
        );
        let mut changed_schema_material = previous
            .expect("fixture previous digest")
            .as_bytes()
            .to_vec();
        changed_schema_material.extend_from_slice(changed_schema_body.as_bytes());
        assert_ne!(
            current,
            fs_blake3::hash_domain(CHECKPOINT_CHAIN_DOMAIN, &changed_schema_material),
            "the checkpoint body schema is semantic"
        );

        let canonical_body = checkpoint_chain_body_json(&input);
        let mut changed_layout_material = previous
            .expect("fixture previous digest")
            .as_bytes()
            .to_vec();
        changed_layout_material.extend_from_slice(
            canonical_body
                .replacen("{\"schema_hex\":", "{ \"schema_hex\":", 1)
                .as_bytes(),
        );
        assert_ne!(
            current,
            fs_blake3::hash_domain(CHECKPOINT_CHAIN_DOMAIN, &changed_layout_material),
            "the canonical JSON layout is semantic"
        );

        let mut changed_framing_material = vec![1];
        changed_framing_material
            .extend_from_slice(previous.expect("fixture previous digest").as_bytes());
        changed_framing_material.extend_from_slice(canonical_body.as_bytes());
        assert_ne!(
            current,
            fs_blake3::hash_domain(CHECKPOINT_CHAIN_DOMAIN, &changed_framing_material),
            "the previous-digest framing is semantic"
        );
    }

    #[test]
    fn staleness_checkpoint_chain_identity_versions_fail_closed() {
        assert_eq!(STALENESS_CHECKPOINT_CHAIN_IDENTITY_VERSION, 2);
        assert!(CHECKPOINT_CHAIN_DOMAIN.ends_with(".v2"));
        assert_eq!(CHECKPOINT_SCHEMA, "fs-roofline-staleness-checkpoint-v2");
        assert_eq!(CHECKPOINT_KERNEL_PREFIX, "roofline-staleness-checkpoint:");
        assert!(CHECKPOINT_SHAPE_PREFIX.contains("ckpt-v2"));
        assert!(LEGACY_CHECKPOINT_SHAPE_PREFIX.contains("ckpt-v1"));
        let entries = vec![checkpoint_identity_entry(1)];
        let input = StalenessCheckpointChainIdentityInput {
            kernel: "kernel-a",
            version: "v1",
            ordinal: 0,
            previous: None,
            entries: &entries,
        };
        assert_ne!(
            staleness_checkpoint_chain_receipt(&input),
            staleness_checkpoint_chain_receipt_with_domain(
                "org.frankensim.fs-roofline.staleness-checkpoint-chain.v3",
                &input,
            )
        );
    }

    #[test]
    fn checkpoint_chain_string_fields_are_injective_and_canonical() {
        let mut entries = vec![checkpoint_identity_entry(1)];
        entries[0].shape_class = "shape\"},\"version\":\"trap".to_string();
        let left = StalenessCheckpointChainIdentityInput {
            kernel: "a\",\"version\":\"b",
            version: "c",
            ordinal: 0,
            previous: None,
            entries: &entries,
        };
        let right = StalenessCheckpointChainIdentityInput {
            kernel: "a",
            version: "b\",\"version\":\"c",
            ..left
        };
        let left_body = checkpoint_chain_body_json(&left);
        let right_body = checkpoint_chain_body_json(&right);
        assert_ne!(
            left_body, right_body,
            "hex framing prevents delimiter collisions"
        );
        let parsed = parse_body(&left_body).expect("canonical hex body round trips");
        assert_eq!(parsed.kernel, left.kernel);
        assert_eq!(parsed.version, left.version);
        assert_eq!(parsed.entries, entries);

        let simple = StalenessCheckpointChainIdentityInput {
            kernel: "j",
            version: "v:1",
            ordinal: 7,
            previous: None,
            entries: &entries,
        };
        let body = checkpoint_chain_body_json(&simple);
        for noncanonical in [
            body.replace("\"kernel_hex\":\"6a\"", "\"kernel_hex\":\"6A\""),
            body.replace("\"kernel_hex\":\"6a\"", "\"kernel_hex\":\"6\""),
            body.replace("\"kernel_hex\":\"6a\"", "\"kernel_hex\":\"6g\""),
        ] {
            assert!(
                parse_body(&noncanonical).is_none(),
                "uppercase, odd, and non-hex transports fail closed"
            );
        }
        assert_eq!(
            checkpoint_ordinal(
                "roofline-ckpt-v2:v:1:00000007",
                CHECKPOINT_SHAPE_PREFIX,
                "v:1"
            ),
            Ok(Some(7))
        );
        assert_eq!(
            checkpoint_ordinal(
                "roofline-ckpt-v2:v:10:00000007",
                CHECKPOINT_SHAPE_PREFIX,
                "v:1"
            ),
            Ok(None),
            "a longer colon-bearing version must not prefix-match"
        );
    }

    #[test]
    fn legacy_valid_rewrite_or_omission_normalizes_to_a_tombstone() {
        let first = checkpoint_identity_entry(1);
        let second = checkpoint_identity_entry(2);

        let mut rewritten = std::collections::BTreeMap::new();
        assert!(normalize_legacy_entries(
            &mut rewritten,
            &[first.clone(), second.clone()]
        ));
        let changed = CheckpointEntry {
            row_hash: fs_blake3::ContentHash([0xA5; 32]),
            ..first.clone()
        };
        assert!(normalize_legacy_entries(
            &mut rewritten,
            &[changed, second.clone()]
        ));
        let rewritten_first = rewritten.get(&first.shape_class).expect("retained rewrite");
        assert_eq!(rewritten_first.verdict, CheckpointVerdict::Corrupt);
        assert_eq!(rewritten_first.row_hash, first.row_hash);
        assert!(corrupt_entry_is_canonical(rewritten_first));

        let mut omitted = std::collections::BTreeMap::new();
        assert!(normalize_legacy_entries(
            &mut omitted,
            &[first.clone(), second.clone()]
        ));
        assert!(normalize_legacy_entries(
            &mut omitted,
            std::slice::from_ref(&second)
        ));
        let omitted_first = omitted.get(&first.shape_class).expect("retained omission");
        assert_eq!(omitted_first.verdict, CheckpointVerdict::Corrupt);
        assert_eq!(omitted_first.row_hash, first.row_hash);
        assert!(corrupt_entry_is_canonical(omitted_first));
    }

    fn temp_db(tag: &str) -> String {
        static NEXT: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
        let n = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        std::env::temp_dir()
            .join(format!(
                "fs-roofline-ckpt-{tag}-{}-{n}.db",
                std::process::id()
            ))
            .display()
            .to_string()
    }

    fn cleanup_db(path: &str) {
        for suffix in ["", "-wal", "-shm", ".fsqlite-wal", ".fsqlite-shm"] {
            let _ = std::fs::remove_file(format!("{path}{suffix}"));
        }
    }

    fn synthetic_axes(fingerprint: u64) -> MachineAxes {
        // Roofs far above any real machine (bead xjhz): cache-resident test
        // kernels must never outrun the fixture roof.
        MachineAxes {
            fingerprint,
            cpu_brand: "synthetic".to_string(),
            logical_cpus: 8,
            bandwidth_single_gbs: 100_000.0,
            bandwidth_all_core_gbs: 400_000.0,
            peak_single_gflops: 50_000.0,
            peak_all_core_gflops: 300_000.0,
        }
    }

    fn trusted_baseline(axes: &MachineAxes) -> (BaselineAxes, BaselineIdentity) {
        let identity =
            BaselineIdentity::current(axes, "test-firmware").expect("valid synthetic identity");
        let candidates: Vec<_> = (0_u64..3)
            .map(|ordinal| {
                BaselineCandidate::from_receipt(
                    axes.clone(),
                    identity.clone(),
                    fs_blake3::hash_domain(
                        "fs-roofline.checkpoint-baseline-source.v1",
                        &ordinal.to_le_bytes(),
                    ),
                )
                .expect("valid synthetic candidate")
            })
            .collect();
        let baseline = promote_baseline(
            &candidates,
            "test-operator",
            "deterministic checkpoint fixture",
            test_day().saturating_sub(10),
            90,
        )
        .expect("valid synthetic baseline");
        (baseline, identity)
    }

    fn attested_policy(
        baseline: &BaselineAxes,
        identity: &BaselineIdentity,
    ) -> AttestedAxisBaselinePolicy {
        AttestedAxisBaselinePolicy::from_verified(
            baseline.clone(),
            identity.clone(),
            test_day(),
            crate::PromotionAttestation::new("test-authority", "test-signature"),
            baseline.provenance().source_receipts().to_vec(),
            crate::PromotionAuthorityDecision::new(
                crate::KeyVerdict::Authorized,
                fs_blake3::hash_domain(
                    "fs-roofline.checkpoint-test-policy.v1",
                    baseline.content_hash().as_bytes(),
                ),
            ),
        )
    }

    fn test_day() -> u64 {
        crate::days_since_epoch_now().expect("unit-test clock after Unix epoch")
    }

    const CONFIG: ProductionRunConfig = ProductionRunConfig {
        n: 1 << 10,
        warmup: 0,
        reps: 1,
    };

    const TEST_DEPGRAPH_RECEIPT: &str = "{\"schema\":\"fs-roofline-synthetic-dependency-receipt-v1\",\"purpose\":\"checkpoint-battery\"}";

    fn test_binding() -> DependencyReceiptBinding {
        let digest = fs_blake3::hash_domain(
            fs_session::GEMM_DEPGRAPH_RECEIPT_DOMAIN,
            TEST_DEPGRAPH_RECEIPT.as_bytes(),
        );
        DependencyReceiptBinding::from_parts(TEST_DEPGRAPH_RECEIPT, digest)
            .expect("test receipt digest agrees")
    }

    struct Fixture {
        ledger: Ledger,
        baseline: BaselineAxes,
        kernels: Vec<(String, String)>,
        recorded_at: i64,
    }

    /// Record one sealed receipt-backed production run into `ledger`.
    fn record_one_run(ledger: &Ledger) -> (Vec<(String, String)>, i64) {
        let axes = synthetic_axes(FINGERPRINT);
        let (baseline, identity) = trusted_baseline(&axes);
        let policy = attested_policy(&baseline, &identity);
        let probe = ProductionProbe::from_observed(axes.clone());
        let post = axes.clone();
        let run = probe
            .run_with_test_receipt(
                CONFIG,
                policy,
                default_registry(1 << 10).expect("bounded registry fixture"),
                move || post,
                TEST_DEPGRAPH_RECEIPT,
            )
            .expect("sealed checkpoint fixture");
        assert!(run.citation_eligible());
        let kernels = run
            .results()
            .iter()
            .map(|result| (result.kernel.clone(), result.version.clone()))
            .collect();
        let op = run.record(ledger).expect("record checkpoint fixture");
        let recorded_at = ledger
            .op(op)
            .unwrap()
            .expect("recorded op")
            .t_end
            .expect("finished op");
        (kernels, recorded_at)
    }

    fn fixture(db: &str) -> Fixture {
        let ledger = Ledger::open(db).expect("open ledger");
        let (kernels, recorded_at) = record_one_run(&ledger);
        let axes = synthetic_axes(FINGERPRINT);
        let (baseline, _) = trusted_baseline(&axes);
        Fixture {
            ledger,
            baseline,
            kernels,
            recorded_at,
        }
    }

    fn seal(fx: &Fixture) -> CheckpointReceipt {
        checkpoint_history_with_dependency(
            &fx.ledger,
            &fx.kernels[0].0,
            &fx.kernels[0].1,
            FINGERPRINT,
            fx.baseline.content_hash(),
            Some(test_binding()),
        )
        .expect("seal checkpoint")
    }

    fn seal_legacy(fx: &Fixture) -> (fs_blake3::ContentHash, Vec<CheckpointEntry>) {
        let (kernel, version) = &fx.kernels[0];
        let matching = match select_matching_rows(
            &fx.ledger,
            kernel,
            version,
            FINGERPRINT,
            Some(fx.baseline.content_hash()),
        )
        .expect("select rows for legacy fixture")
        {
            RowSelection::Rows(rows) => rows,
            RowSelection::Verdict(verdict) => {
                panic!("legacy fixture unexpectedly classified {verdict:?}")
            }
        };
        let mut entries = matching
            .iter()
            .map(|row| {
                let validated = validate_roofline_row(
                    &fx.ledger,
                    row,
                    kernel,
                    version,
                    FINGERPRINT,
                    fx.baseline.content_hash(),
                    Some(test_binding()),
                )
                .expect("validate legacy fixture row");
                seal_entry(row, validated)
            })
            .collect::<Vec<_>>();
        entries.sort_by(|left, right| left.shape_class.cmp(&right.shape_class));
        let body = legacy_body_json(kernel, version, 0, None, &entries);
        let digest = legacy_chain_digest(None, &body);
        fx.ledger
            .tune_put_if_absent(
                &checkpoint_kernel(kernel),
                &legacy_checkpoint_shape_class(version, 0),
                &roofline_machine_key(FINGERPRINT, fx.baseline.content_hash()),
                &legacy_checkpoint_params(0, digest),
                &body,
            )
            .expect("insert valid legacy checkpoint");
        (digest, entries)
    }

    fn fast(fx: &Fixture, kernel: &str, version: &str, at: i64) -> Staleness {
        staleness_at_checkpointed_with(
            &fx.ledger,
            kernel,
            version,
            FINGERPRINT,
            Some(fx.baseline.content_hash()),
            at,
            crate::executable_build_identity().expect("build identity"),
            Some(test_binding()),
        )
        .expect("checkpointed probe")
    }

    fn exhaustive(fx: &Fixture, kernel: &str, version: &str, at: i64) -> Staleness {
        crate::staleness_at_with_build_and_dependency(
            &fx.ledger,
            kernel,
            version,
            FINGERPRINT,
            Some(fx.baseline.content_hash()),
            at,
            crate::executable_build_identity().expect("build identity"),
            Some(test_binding()),
        )
        .expect("exhaustive probe")
    }

    /// The stored roofline row for one production kernel.
    fn production_row(ledger: &Ledger, kernel: &str) -> fs_ledger::TuneRow {
        let mut rows: Vec<_> = ledger
            .tune_rows(kernel)
            .expect("tune rows")
            .into_iter()
            .filter(|row| row.shape_class.contains(":run="))
            .collect();
        assert!(!rows.is_empty(), "expected roofline rows for {kernel}");
        rows.pop().expect("row")
    }

    #[test]
    fn checkpointed_verdicts_match_exhaustive_across_the_lattice() {
        let db = temp_db("equivalence");
        let fx = fixture(&db);
        let receipt = seal(&fx);
        assert_eq!(receipt.ordinal, 0);
        assert_eq!(receipt.corrupt_rows, 0);
        assert!(receipt.rows >= 1);

        let (kernel, version) = &fx.kernels[0];
        // Fresh / Expired / ClockRollback classification parity.
        for at in [
            fx.recorded_at + 1,
            fx.recorded_at + STALENESS_MAX_AGE_NS,
            fx.recorded_at + STALENESS_MAX_AGE_NS + 1,
            fx.recorded_at - 1,
        ] {
            assert_eq!(
                fast(&fx, kernel, version, at),
                exhaustive(&fx, kernel, version, at),
                "fast/exhaustive divergence at offset {}",
                at - fx.recorded_at
            );
        }
        assert_eq!(
            fast(&fx, kernel, version, fx.recorded_at + 1),
            Staleness::Fresh
        );
        assert_eq!(
            fast(
                &fx,
                kernel,
                version,
                fx.recorded_at + STALENESS_MAX_AGE_NS + 1
            ),
            Staleness::Expired
        );
        assert_eq!(
            fast(&fx, kernel, version, fx.recorded_at - 1),
            Staleness::ClockRollback
        );

        // Lattice-prefix verdicts decided before the chain is even loaded.
        assert_eq!(
            fast(&fx, "never-measured-kernel", version, fx.recorded_at + 1),
            Staleness::NeverMeasured
        );
        let wrong_fp = staleness_at_checkpointed_with(
            &fx.ledger,
            kernel,
            version,
            FINGERPRINT + 1,
            Some(fx.baseline.content_hash()),
            fx.recorded_at + 1,
            crate::executable_build_identity().expect("build identity"),
            Some(test_binding()),
        )
        .expect("wrong-fingerprint probe");
        assert_eq!(wrong_fp, Staleness::FingerprintDrift);
        let no_baseline = staleness_at_checkpointed_with(
            &fx.ledger,
            kernel,
            version,
            FINGERPRINT,
            None,
            fx.recorded_at + 1,
            crate::executable_build_identity().expect("build identity"),
            Some(test_binding()),
        )
        .expect("no-baseline probe");
        assert_eq!(no_baseline, Staleness::BaselineUnavailable);

        // Foreign build (injected): both paths must agree on BuildDrift.
        let foreign_build = fs_blake3::hash_domain("fs-roofline.ckpt-test-build.v1", b"other");
        let fast_foreign = staleness_at_checkpointed_with(
            &fx.ledger,
            kernel,
            version,
            FINGERPRINT,
            Some(fx.baseline.content_hash()),
            fx.recorded_at + 1,
            foreign_build,
            Some(test_binding()),
        )
        .expect("foreign-build fast probe");
        let exhaustive_foreign = crate::staleness_at_with_build_and_dependency(
            &fx.ledger,
            kernel,
            version,
            FINGERPRINT,
            Some(fx.baseline.content_hash()),
            fx.recorded_at + 1,
            foreign_build,
            Some(test_binding()),
        )
        .expect("foreign-build exhaustive probe");
        assert_eq!(fast_foreign, exhaustive_foreign);
        assert_eq!(fast_foreign, Staleness::BuildDrift);
        cleanup_db(&db);
    }

    #[test]
    fn malformed_legacy_v1_checkpoint_falls_back_to_exhaustive() {
        let db = temp_db("legacy-v1-malformed");
        let fx = fixture(&db);
        let (kernel, version) = &fx.kernels[0];
        fx.ledger
            .tune_put(
                &format!("roofline-staleness-checkpoint:{kernel}"),
                &format!("roofline-ckpt-v1:{version}:0"),
                &roofline_machine_key(FINGERPRINT, fx.baseline.content_hash()),
                "{}",
                "{}",
            )
            .expect("insert legacy checkpoint row");

        let at = fx.recorded_at + 1;
        let verdict = fast(&fx, kernel, version, at);
        assert_eq!(
            verdict,
            exhaustive(&fx, kernel, version, at),
            "an unverifiable v1 checkpoint must never authorize the fast path"
        );
        assert_eq!(verdict, Staleness::Fresh);
        cleanup_db(&db);
    }

    #[test]
    fn legacy_v1_chain_is_bound_into_v2_genesis() {
        let db = temp_db("legacy-v1-bridge");
        let fx = fixture(&db);
        let (kernel, version) = &fx.kernels[0];
        let (legacy_tip, legacy_entries) = seal_legacy(&fx);
        assert_eq!(
            fast(&fx, kernel, version, fx.recorded_at + 1),
            Staleness::Fresh,
            "a strictly verified legacy seal remains active before migration"
        );

        let receipt = seal(&fx);
        assert_eq!(receipt.ordinal, 0, "v2 ordinals restart at zero");
        let current = fx
            .ledger
            .tune_get(
                &checkpoint_kernel(kernel),
                &checkpoint_shape_class(version, 0),
                &roofline_machine_key(FINGERPRINT, fx.baseline.content_hash()),
            )
            .expect("read v2 genesis")
            .expect("v2 genesis exists");
        let body = parse_body(&current.measured).expect("v2 genesis parses");
        assert_eq!(body.prev, Some(legacy_tip));
        assert_eq!(body.entries, legacy_entries);
        assert_eq!(
            fast(&fx, kernel, version, fx.recorded_at + 1),
            Staleness::Fresh
        );
        cleanup_db(&db);
    }

    #[test]
    fn historical_v1_tombstone_reseals_migrate_monotonically() {
        let db = temp_db("legacy-v1-tombstone-migration");
        let fx = fixture(&db);
        let (kernel, version) = &fx.kernels[0];
        let (valid_tip, valid_entries) = seal_legacy(&fx);
        let original = production_row(&fx.ledger, kernel);
        let original_entry = valid_entries
            .iter()
            .find(|entry| entry.shape_class == original.shape_class)
            .expect("v1 genesis covers the production row")
            .clone();

        let forged = original
            .measured
            .replace("\"dispersion\":", "\"dispersion\": ");
        assert_ne!(forged, original.measured);
        fx.ledger
            .tune_put(
                &original.kernel,
                &original.shape_class,
                &original.machine,
                &original.params,
                &forged,
            )
            .expect("tamper retained production row");
        let tampered_row = production_row(&fx.ledger, kernel);
        let mut corrupt_entries = valid_entries.clone();
        let corrupt_index = corrupt_entries
            .iter()
            .position(|entry| entry.shape_class == tampered_row.shape_class)
            .expect("tampered row remains covered");
        corrupt_entries[corrupt_index] = seal_entry(&tampered_row, None);
        assert_ne!(
            corrupt_entries[corrupt_index].row_hash, original_entry.row_hash,
            "historical v1 recomputed the tombstone hash from tampered bytes"
        );
        let corrupt_body = legacy_body_json(kernel, version, 1, Some(valid_tip), &corrupt_entries);
        let corrupt_tip = legacy_chain_digest(Some(valid_tip), &corrupt_body);
        fx.ledger
            .tune_put_if_absent(
                &checkpoint_kernel(kernel),
                &legacy_checkpoint_shape_class(version, 1),
                &roofline_machine_key(FINGERPRINT, fx.baseline.content_hash()),
                &legacy_checkpoint_params(1, corrupt_tip),
                &corrupt_body,
            )
            .expect("insert historical v1 corrupt snapshot");

        fx.ledger
            .tune_put(
                &original.kernel,
                &original.shape_class,
                &original.machine,
                &original.params,
                &original.measured,
            )
            .expect("restore retained production row");
        let restored_row = production_row(&fx.ledger, kernel);
        let mut restored_entries = corrupt_entries.clone();
        restored_entries[corrupt_index] = seal_entry(&restored_row, None);
        assert_ne!(
            restored_entries[corrupt_index].row_hash, corrupt_entries[corrupt_index].row_hash,
            "historical v1 recomputed an inherited tombstone after restoration"
        );
        let restored_body =
            legacy_body_json(kernel, version, 2, Some(corrupt_tip), &restored_entries);
        let restored_tip = legacy_chain_digest(Some(corrupt_tip), &restored_body);
        fx.ledger
            .tune_put_if_absent(
                &checkpoint_kernel(kernel),
                &legacy_checkpoint_shape_class(version, 2),
                &roofline_machine_key(FINGERPRINT, fx.baseline.content_hash()),
                &legacy_checkpoint_params(2, restored_tip),
                &restored_body,
            )
            .expect("insert historical v1 restored-row tombstone snapshot");

        let migrated = seal(&fx);
        assert_eq!(migrated.ordinal, 0, "v2 migration starts at ordinal zero");
        assert!(migrated.corrupt_rows >= 1);
        let current = fx
            .ledger
            .tune_get(
                &checkpoint_kernel(kernel),
                &checkpoint_shape_class(version, 0),
                &roofline_machine_key(FINGERPRINT, fx.baseline.content_hash()),
            )
            .expect("read v2 migration row")
            .expect("v2 migration row exists");
        let body = parse_body(&current.measured).expect("v2 migration row parses");
        let migrated_entry = body
            .entries
            .iter()
            .find(|entry| entry.shape_class == original.shape_class)
            .expect("v2 migration retains the historical row");
        assert_eq!(migrated_entry.verdict, CheckpointVerdict::Corrupt);
        assert_eq!(
            migrated_entry.row_hash, original_entry.row_hash,
            "migration pins the first incident to the pre-tamper covered identity"
        );
        assert_eq!(
            fast(&fx, kernel, version, fx.recorded_at + 1),
            Staleness::CorruptEvidence,
            "v1 restoration cannot erase a tombstone during migration"
        );
        cleanup_db(&db);
    }

    #[test]
    fn late_legacy_append_breaks_the_v2_bridge() {
        let db = temp_db("legacy-v1-late-append");
        let fx = fixture(&db);
        let (kernel, version) = &fx.kernels[0];
        let (legacy_tip, entries) = seal_legacy(&fx);
        seal(&fx);

        let late_body = legacy_body_json(kernel, version, 1, Some(legacy_tip), &entries);
        let late_digest = legacy_chain_digest(Some(legacy_tip), &late_body);
        fx.ledger
            .tune_put_if_absent(
                &checkpoint_kernel(kernel),
                &legacy_checkpoint_shape_class(version, 1),
                &roofline_machine_key(FINGERPRINT, fx.baseline.content_hash()),
                &legacy_checkpoint_params(1, late_digest),
                &late_body,
            )
            .expect("append late legacy checkpoint");

        assert_eq!(
            fast(&fx, kernel, version, fx.recorded_at + 1),
            exhaustive(&fx, kernel, version, fx.recorded_at + 1),
            "a moved legacy tip breaks v2 continuity and forces exhaustive validation"
        );
        let refusal = checkpoint_history_with_dependency(
            &fx.ledger,
            kernel,
            version,
            FINGERPRINT,
            fx.baseline.content_hash(),
            Some(test_binding()),
        )
        .expect_err("a late legacy append must block v2 extension");
        assert!(matches!(refusal, LedgerError::Invalid { .. }));
        cleanup_db(&db);
    }

    #[test]
    fn fast_path_costs_two_reads_and_undercuts_exhaustive() {
        let db = temp_db("budget");
        let fx = fixture(&db);
        seal(&fx);
        let (kernel, version) = &fx.kernels[0];

        let before = fx.ledger.read_queries();
        assert_eq!(
            fast(&fx, kernel, version, fx.recorded_at + 1),
            Staleness::Fresh
        );
        let fast_reads = fx.ledger.read_queries() - before;

        let before = fx.ledger.read_queries();
        assert_eq!(
            exhaustive(&fx, kernel, version, fx.recorded_at + 1),
            Staleness::Fresh
        );
        let exhaustive_reads = fx.ledger.read_queries() - before;

        // Fully covered fast path: one tune_rows for the production kernel,
        // one for the checkpoint chain. Nothing per-row.
        assert_eq!(
            fast_reads, 2,
            "covered fast path must not scale with history"
        );
        assert!(
            exhaustive_reads > fast_reads,
            "exhaustive ({exhaustive_reads}) must cost more than checkpointed ({fast_reads})"
        );
        cleanup_db(&db);
    }

    #[test]
    fn tampered_production_row_is_corrupt_under_the_checkpoint() {
        let db = temp_db("tamper");
        let fx = fixture(&db);
        seal(&fx);
        let (kernel, version) = &fx.kernels[0];
        assert_eq!(
            fast(&fx, kernel, version, fx.recorded_at + 1),
            Staleness::Fresh
        );

        let row = production_row(&fx.ledger, kernel);
        let forged = row.measured.replace("\"dispersion\":", "\"dispersion\": ");
        assert_ne!(forged, row.measured);
        fx.ledger
            .tune_put(
                &row.kernel,
                &row.shape_class,
                &row.machine,
                &row.params,
                &forged,
            )
            .expect("overwrite row");

        assert_eq!(
            fast(&fx, kernel, version, fx.recorded_at + 1),
            Staleness::CorruptEvidence,
            "content-hash mismatch against the sealed entry must be corrupt"
        );
        let checkpoint_rows_before = fx
            .ledger
            .tune_rows(&checkpoint_kernel(kernel))
            .expect("checkpoint rows before refused reseal")
            .len();
        let reseal = checkpoint_history_with_dependency(
            &fx.ledger,
            kernel,
            version,
            FINGERPRINT,
            fx.baseline.content_hash(),
            Some(test_binding()),
        )
        .expect_err("resealing must not bless changed previously valid bytes");
        assert!(matches!(reseal, LedgerError::Invalid { .. }));
        assert!(
            !fx.ledger.in_transaction(),
            "owned failure must leave no open transaction"
        );
        assert_eq!(
            fx.ledger
                .tune_rows(&checkpoint_kernel(kernel))
                .expect("checkpoint rows after refused reseal")
                .len(),
            checkpoint_rows_before,
            "owned rollback must leave the checkpoint chain unchanged"
        );
        cleanup_db(&db);
    }

    #[test]
    fn checkpoint_sealing_refuses_a_caller_owned_transaction() {
        let db = temp_db("caller-owned-transaction");
        let fx = fixture(&db);
        let (kernel, version) = &fx.kernels[0];
        fx.ledger.begin().expect("begin caller transaction");

        let error = checkpoint_history_with_dependency(
            &fx.ledger,
            kernel,
            version,
            FINGERPRINT,
            fx.baseline.content_hash(),
            Some(test_binding()),
        )
        .expect_err("checkpoint sealing must own its rollback boundary");
        assert!(matches!(
            &error,
            LedgerError::Invalid { field, problem }
                if field == "staleness-checkpoint-transaction"
                    && problem.contains("caller transaction")
        ));
        assert!(
            fx.ledger.in_transaction(),
            "the refused call must not commit or roll back caller state"
        );
        assert!(
            fx.ledger
                .tune_rows(&checkpoint_kernel(kernel))
                .expect("read checkpoint rows inside caller transaction")
                .is_empty(),
            "refusal happens before any checkpoint insert"
        );
        fx.ledger.rollback().expect("roll back caller transaction");
        cleanup_db(&db);
    }

    #[test]
    fn covered_row_missing_from_the_ledger_is_corrupt() {
        // Rollback simulation: rebuild a ledger holding the checkpoint chain
        // and the second run's rows but MISSING the first run's row — as if
        // history was truncated underneath the seal.
        let db = temp_db("rollback-src");
        let fx = fixture(&db);
        let (kernels2, recorded_at2) = record_one_run(&fx.ledger);
        assert_eq!(fx.kernels, kernels2);
        seal(&fx); // covers both runs' rows
        let (kernel, version) = &fx.kernels[0];
        assert_eq!(
            fast(&fx, kernel, version, recorded_at2 + 1),
            Staleness::Fresh
        );

        let replay_db = temp_db("rollback-dst");
        let replay = Ledger::open(&replay_db).expect("open replay ledger");
        // Copy the checkpoint chain verbatim.
        for row in fx
            .ledger
            .tune_rows(&checkpoint_kernel(kernel))
            .expect("chain rows")
        {
            replay
                .tune_put(
                    &row.kernel,
                    &row.shape_class,
                    &row.machine,
                    &row.params,
                    &row.measured,
                )
                .expect("copy chain row");
        }
        // Copy only the NEWER of the two production rows.
        let mut rows: Vec<_> = fx
            .ledger
            .tune_rows(kernel)
            .expect("production rows")
            .into_iter()
            .filter(|row| row.shape_class.contains(":run="))
            .collect();
        assert_eq!(rows.len(), 2, "two runs leave two rows");
        rows.sort_by_key(|row| row.shape_class.clone());
        let kept = rows.pop().expect("newest row");
        replay
            .tune_put(
                &kept.kernel,
                &kept.shape_class,
                &kept.machine,
                &kept.params,
                &kept.measured,
            )
            .expect("copy surviving row");

        let verdict = staleness_at_checkpointed_with(
            &replay,
            kernel,
            version,
            FINGERPRINT,
            Some(fx.baseline.content_hash()),
            recorded_at2 + 1,
            crate::executable_build_identity().expect("build identity"),
            Some(test_binding()),
        )
        .expect("rollback probe");
        assert_eq!(
            verdict,
            Staleness::CorruptEvidence,
            "a sealed row vanishing from history must be corrupt, not silently fresh"
        );
        let reseal = checkpoint_history_with_dependency(
            &replay,
            kernel,
            version,
            FINGERPRINT,
            fx.baseline.content_hash(),
            Some(test_binding()),
        )
        .expect_err("resealing must not erase a previously valid missing row");
        assert!(matches!(reseal, LedgerError::Invalid { .. }));
        cleanup_db(&db);
        cleanup_db(&replay_db);
    }

    #[test]
    fn removal_of_every_covered_production_row_is_corrupt() {
        let db = temp_db("total-rollback-src");
        let fx = fixture(&db);
        seal(&fx);
        let (kernel, version) = &fx.kernels[0];

        let replay_db = temp_db("total-rollback-dst");
        let replay = Ledger::open(&replay_db).expect("open total-rollback ledger");
        for row in fx
            .ledger
            .tune_rows(&checkpoint_kernel(kernel))
            .expect("checkpoint rows")
        {
            replay
                .tune_put(
                    &row.kernel,
                    &row.shape_class,
                    &row.machine,
                    &row.params,
                    &row.measured,
                )
                .expect("copy checkpoint row");
        }

        let verdict = staleness_at_checkpointed_with(
            &replay,
            kernel,
            version,
            FINGERPRINT,
            Some(fx.baseline.content_hash()),
            fx.recorded_at + 1,
            crate::executable_build_identity().expect("build identity"),
            Some(test_binding()),
        )
        .expect("total rollback probe");
        assert_eq!(
            verdict,
            Staleness::CorruptEvidence,
            "verified nonempty history must be consulted before NeverMeasured"
        );
        cleanup_db(&db);
        cleanup_db(&replay_db);
    }

    #[test]
    fn tampered_chain_fails_closed_to_exhaustive_and_blocks_sealing() {
        let db = temp_db("chain-tamper");
        let fx = fixture(&db);
        seal(&fx);
        let (kernel, version) = &fx.kernels[0];

        // Flip bytes inside the sealed body. The chain no longer verifies.
        let chain_kernel = checkpoint_kernel(kernel);
        let row = fx
            .ledger
            .tune_rows(&chain_kernel)
            .expect("chain rows")
            .pop()
            .expect("chain row");
        let forged = row
            .measured
            .replace("\"verdict\":\"valid\"", "\"verdict\":\"corrupt\"");
        assert_ne!(
            forged, row.measured,
            "fixture must have a valid entry to flip"
        );
        fx.ledger
            .tune_put(
                &row.kernel,
                &row.shape_class,
                &row.machine,
                &row.params,
                &forged,
            )
            .expect("tamper chain");

        // Fail closed: the fast path falls back to the exhaustive verdict.
        let before = fx.ledger.read_queries();
        assert_eq!(
            fast(&fx, kernel, version, fx.recorded_at + 1),
            Staleness::Fresh
        );
        let fallback_reads = fx.ledger.read_queries() - before;
        assert!(
            fallback_reads > 2,
            "broken chain must trigger the exhaustive fallback, not the covered path"
        );

        // And extending a broken chain is refused outright.
        let refusal = checkpoint_history_with_dependency(
            &fx.ledger,
            kernel,
            version,
            FINGERPRINT,
            fx.baseline.content_hash(),
            Some(test_binding()),
        )
        .expect_err("sealing over a broken chain must refuse");
        assert!(matches!(refusal, LedgerError::Invalid { .. }));
        cleanup_db(&db);
    }

    #[test]
    fn ordinal_gaps_and_duplicate_ordinals_break_the_chain() {
        let db = temp_db("gap");
        let fx = fixture(&db);
        seal(&fx);
        let (kernel, version) = &fx.kernels[0];
        let chain_kernel = checkpoint_kernel(kernel);
        let row = fx
            .ledger
            .tune_rows(&chain_kernel)
            .expect("chain rows")
            .pop()
            .expect("chain row");

        // Insert a structurally plausible row at a gapped ordinal (2 with no 1).
        fx.ledger
            .tune_put(
                &row.kernel,
                &checkpoint_shape_class(version, 2),
                &row.machine,
                &row.params,
                &row.measured,
            )
            .expect("gapped row");
        let before = fx.ledger.read_queries();
        assert_eq!(
            fast(&fx, kernel, version, fx.recorded_at + 1),
            Staleness::Fresh
        );
        assert!(
            fx.ledger.read_queries() - before > 2,
            "gapped chain must fail closed to exhaustive"
        );
        cleanup_db(&db);
    }

    #[test]
    fn tombstones_survive_row_restoration_and_resealing() {
        let db = temp_db("tombstone");
        let fx = fixture(&db);
        let (kernel, version) = &fx.kernels[0];
        let original = production_row(&fx.ledger, kernel);

        // Tamper, then seal: the exhaustive pass records tombstones.
        let forged = original
            .measured
            .replace("\"dispersion\":", "\"dispersion\": ");
        assert_ne!(forged, original.measured);
        fx.ledger
            .tune_put(
                &original.kernel,
                &original.shape_class,
                &original.machine,
                &original.params,
                &forged,
            )
            .expect("tamper row");
        let sealed = seal(&fx);
        assert!(
            sealed.corrupt_rows >= 1,
            "tampered history must seal tombstones"
        );

        // Restore the original bytes: the row would re-verify exhaustively...
        fx.ledger
            .tune_put(
                &original.kernel,
                &original.shape_class,
                &original.machine,
                &original.params,
                &original.measured,
            )
            .expect("restore row");
        assert_eq!(
            exhaustive(&fx, kernel, version, fx.recorded_at + 1),
            Staleness::Fresh,
            "precondition: restored bytes verify exhaustively"
        );
        // ...but the tombstone is permanent on the fast path...
        assert_eq!(
            fast(&fx, kernel, version, fx.recorded_at + 1),
            Staleness::CorruptEvidence,
            "restoration must never un-corrupt a sealed tombstone"
        );
        // ...and re-sealing inherits it rather than forgiving it.
        let resealed = seal(&fx);
        assert_eq!(resealed.ordinal, sealed.ordinal + 1);
        assert!(
            resealed.corrupt_rows >= sealed.corrupt_rows,
            "re-checkpointing must preserve every prior tombstone"
        );
        assert_eq!(
            fast(&fx, kernel, version, fx.recorded_at + 1),
            Staleness::CorruptEvidence
        );
        cleanup_db(&db);
    }

    #[test]
    fn delta_rows_recorded_after_the_seal_still_classify() {
        let db = temp_db("delta");
        let fx = fixture(&db);
        seal(&fx);
        let (kernels2, recorded_at2) = record_one_run(&fx.ledger);
        assert_eq!(fx.kernels, kernels2);
        let (kernel, version) = &fx.kernels[0];

        // The second run's row is not covered by the checkpoint; the fast
        // path validates it exhaustively and classifies over BOTH rows.
        assert_eq!(
            fast(&fx, kernel, version, recorded_at2 + 1),
            exhaustive(&fx, kernel, version, recorded_at2 + 1)
        );
        assert_eq!(
            fast(&fx, kernel, version, recorded_at2 + 1),
            Staleness::Fresh
        );

        // Delta validation costs more than covered but the next seal
        // re-covers everything.
        let before = fx.ledger.read_queries();
        let _ = fast(&fx, kernel, version, recorded_at2 + 1);
        let delta_reads = fx.ledger.read_queries() - before;
        assert!(delta_reads > 2, "delta rows must be exhaustively validated");

        let receipt = seal(&fx);
        assert_eq!(receipt.ordinal, 1);
        let before = fx.ledger.read_queries();
        assert_eq!(
            fast(&fx, kernel, version, recorded_at2 + 1),
            Staleness::Fresh
        );
        assert_eq!(
            fx.ledger.read_queries() - before,
            2,
            "after re-sealing, the fast path is fully covered again"
        );
        cleanup_db(&db);
    }

    #[test]
    #[ignore = "scale fixture: ~40 production runs; run explicitly"]
    fn fast_path_read_cost_stays_constant_as_history_grows() {
        let db = temp_db("scale");
        let fx = fixture(&db);
        let (kernel, version) = &fx.kernels[0];
        let mut newest = fx.recorded_at;
        let mut exhaustive_read_counts = Vec::new();
        for _ in 0..39 {
            let (_, at) = record_one_run(&fx.ledger);
            newest = at;
        }
        seal(&fx);

        let before = fx.ledger.read_queries();
        assert_eq!(fast(&fx, kernel, version, newest + 1), Staleness::Fresh);
        let fast_reads = fx.ledger.read_queries() - before;
        assert_eq!(fast_reads, 2, "40 covered rows must still cost two reads");

        let before = fx.ledger.read_queries();
        assert_eq!(
            exhaustive(&fx, kernel, version, newest + 1),
            Staleness::Fresh
        );
        exhaustive_read_counts.push(fx.ledger.read_queries() - before);
        assert!(
            exhaustive_read_counts[0] > 40,
            "exhaustive cost must scale with history ({} reads)",
            exhaustive_read_counts[0]
        );
        cleanup_db(&db);
    }
}

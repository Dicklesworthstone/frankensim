//! Crash-stable mutation claims and terminal receipts for fs-session.
//!
//! A claim is committed before caller-owned work begins. Exact retries then
//! observe one of three durable states: newly claimed authority, an existing
//! indeterminate Pending claim, or a verified terminal receipt. Pending work
//! is never silently executed again.
//!
//! Terminal receipt bytes and their owned global audit events commit in one
//! transaction. Every variable-size value is bounded before materialization,
//! payload and receipt BLOB hashes are computed by this crate, and terminal
//! reads rejoin and rehash every owned event. The registry therefore closes
//! both the database-commit/in-memory-cursor window and ordinary partial-row
//! tampering without treating audit JSON as executable recovery state.

use std::collections::BTreeSet;

use fsqlite::PreparedStatement;

use super::*;

/// Version of the canonical claim, terminal, event-link, and batch envelopes.
pub const SESSION_REGISTRY_ROW_SCHEMA_VERSION: i64 = 1;

/// Maximum UTF-8 bytes in a terminal mutation kind.
pub const MAX_SESSION_TERMINAL_KIND_BYTES: usize = 64;

/// Maximum UTF-8 bytes in an immutable ledger scope.
pub const MAX_SESSION_TERMINAL_SCOPE_BYTES: usize = 128;

/// Maximum bytes in one canonical mutation payload BLOB.
pub const MAX_SESSION_CLAIM_PAYLOAD_BYTES: usize = 1024 * 1024;

/// Maximum bytes in one canonical terminal receipt BLOB.
pub const MAX_SESSION_TERMINAL_RECEIPT_BYTES: usize = 1024 * 1024;

/// Maximum UTF-8 bytes in one owned audit-event kind.
pub const MAX_SESSION_TERMINAL_EVENT_KIND_BYTES: usize = 256;

/// Maximum bytes in one owned audit-event JSON payload.
pub const MAX_SESSION_TERMINAL_EVENT_PAYLOAD_BYTES: usize = 1024 * 1024;

/// Maximum terminal identities admitted by one atomic batch.
pub const MAX_SESSION_FLUSH_TERMINALS: usize = 1024;

/// Maximum canonical batch witnesses retained for one terminal identity.
pub const MAX_SESSION_TERMINAL_BATCH_MEMBERSHIPS: usize = 1024;

/// Maximum claims inspected by one generic generation recovery probe.
pub const MAX_SESSION_RECOVERY_PROBE_CLAIMS: usize = 8192;

/// Maximum submission claims inspected by one pause-generation fence.
pub const MAX_SESSION_PAUSE_FENCE_SUBMISSIONS: usize = 4096;

/// Maximum audit events admitted by one atomic batch.
pub const MAX_SESSION_FLUSH_EVENTS: usize = 1024;

/// Maximum conservatively encoded bytes admitted by one terminal batch.
pub const MAX_SESSION_FLUSH_ENCODED_BYTES: usize = 4 * 1024 * 1024;

const SESSION_CLAIM_ROW_FRAMING_BYTES: usize = 256;
const SESSION_TERMINAL_ROW_FRAMING_BYTES: usize = 96;
const SESSION_EVENT_ROW_FRAMING_BYTES: usize = 64;
const SESSION_CLAIM_HASH_DOMAIN: &[u8] = b"org.frankensim.fs-ledger.session-mutation-claim.v1\0";
// The tracked v6 schema shipped without a wired registry writer. These v2
// domains therefore define the first supported batch/event preimages; rows
// from an earlier uncommitted scaffold are intentionally not auto-trusted.
const SESSION_BATCH_HASH_DOMAIN: &[u8] = b"org.frankensim.fs-ledger.session-flush-batch.v2\0";
const SESSION_EVENTS_HASH_DOMAIN: &[u8] = b"org.frankensim.fs-ledger.session-terminal-events.v2\0";
const PRECLAIM_REQUIRED_SUBMISSION_KIND: &str = "submission";
const PAUSE_ACKNOWLEDGEMENT_KIND: &str = "pause-acknowledgement";

/// Semantic version of a ledger-minted solver-checkpoint receipt.
pub const SOLVER_CHECKPOINT_RECEIPT_IDENTITY_VERSION: u32 = 1;
/// Exact artifact kind required for attested solver-state snapshots.
pub const SOLVER_STATE_ARTIFACT_KIND: &str = "solver-state";
/// Maximum solver-state artifact bytes independently materialized at mint and
/// verification time.
pub const MAX_SOLVER_CHECKPOINT_ARTIFACT_BYTES: u64 = 64 * 1024 * 1024;
const SOLVER_CHECKPOINT_RECEIPT_IDENTITY_DOMAIN: &[u8] =
    b"org.frankensim.fs-ledger.solver-checkpoint-receipt.v1\0";
const SOLVER_CHECKPOINT_RECEIPT_TRANSPORT_BYTES: usize =
    4 + 16 + 8 + 8 + 32 + 8 + 32 + 32 + 8 + 8 + 32;

/// Inputs to one immutable solver-checkpoint receipt.
///
/// The drain report is executor-minted with private fields. The ledger derives
/// the run from it, validates the referenced snapshot envelope independently,
/// and never accepts caller-authored `drained=true` evidence.
#[derive(Debug, Clone, Copy)]
pub struct SolverCheckpointClaim<'a> {
    /// Session whose old execution generation drained.
    pub session: u64,
    /// Domain-separated authority of the exact pending pause request.
    pub pause_authority: ContentHash,
    /// Cancellation-gate generation that was drained.
    pub gate_generation: u64,
    /// Existing content-addressed `solver-state` artifact.
    pub solver_state_artifact: ContentHash,
    /// Executor-originated report minted after all registered workers joined.
    pub drain_report: &'a fs_exec::DrainFinalizeReport,
}

/// Immutable private-field receipt minted by one physical ledger.
///
/// Transport decoding reconstructs only an untrusted typed candidate. A
/// consumer MUST call [`Ledger::verify_solver_checkpoint_receipt`] before
/// treating it as completion evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SolverCheckpointReceipt {
    ledger_instance_id: LedgerInstanceId,
    session: u64,
    run: fs_exec::RunId,
    pause_authority: ContentHash,
    gate_generation: u64,
    solver_state_artifact: ContentHash,
    drain_report_hash: ContentHash,
    registered_workers: u64,
    drained_workers: u64,
    content_hash: ContentHash,
}

impl SolverCheckpointReceipt {
    /// Physical ledger that minted and stores this receipt.
    #[must_use]
    pub const fn ledger_instance_id(self) -> LedgerInstanceId {
        self.ledger_instance_id
    }

    /// Session bound into the receipt.
    #[must_use]
    pub const fn session(self) -> u64 {
        self.session
    }

    /// Logical executor run recovered from the validated snapshot/report pair.
    #[must_use]
    pub const fn run(self) -> fs_exec::RunId {
        self.run
    }

    /// Exact pending pause authority.
    #[must_use]
    pub const fn pause_authority(self) -> ContentHash {
        self.pause_authority
    }

    /// Drained cancellation-gate generation.
    #[must_use]
    pub const fn gate_generation(self) -> u64 {
        self.gate_generation
    }

    /// Content hash of the validated solver-state artifact.
    #[must_use]
    pub const fn solver_state_artifact(self) -> ContentHash {
        self.solver_state_artifact
    }

    /// Executor drain/finalize report identity.
    #[must_use]
    pub const fn drain_report_hash(self) -> ContentHash {
        self.drain_report_hash
    }

    /// Workers admitted into the executor drain boundary.
    #[must_use]
    pub const fn registered_workers(self) -> u64 {
        self.registered_workers
    }

    /// Workers joined before finalization.
    #[must_use]
    pub const fn drained_workers(self) -> u64 {
        self.drained_workers
    }

    /// Domain-separated identity of every semantic receipt field.
    #[must_use]
    pub const fn content_hash(self) -> ContentHash {
        self.content_hash
    }

    /// Fixed-width transport for response-loss recovery and process handoff.
    #[must_use]
    pub fn to_bytes(self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(SOLVER_CHECKPOINT_RECEIPT_TRANSPORT_BYTES);
        bytes.extend_from_slice(&SOLVER_CHECKPOINT_RECEIPT_IDENTITY_VERSION.to_le_bytes());
        bytes.extend_from_slice(&self.ledger_instance_id.as_bytes());
        bytes.extend_from_slice(&self.session.to_le_bytes());
        bytes.extend_from_slice(&self.run.0.to_le_bytes());
        bytes.extend_from_slice(self.pause_authority.as_bytes());
        bytes.extend_from_slice(&self.gate_generation.to_le_bytes());
        bytes.extend_from_slice(self.solver_state_artifact.as_bytes());
        bytes.extend_from_slice(self.drain_report_hash.as_bytes());
        bytes.extend_from_slice(&self.registered_workers.to_le_bytes());
        bytes.extend_from_slice(&self.drained_workers.to_le_bytes());
        bytes.extend_from_slice(self.content_hash.as_bytes());
        bytes
    }

    /// Decode an untrusted fixed-width transport candidate.
    ///
    /// This checks self-consistency only. Ledger membership and artifact
    /// existence are re-earned by [`Ledger::verify_solver_checkpoint_receipt`].
    ///
    /// # Errors
    /// Wrong length/version, malformed drain counts, or identity mismatch.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, LedgerError> {
        if bytes.len() != SOLVER_CHECKPOINT_RECEIPT_TRANSPORT_BYTES {
            return Err(invalid(
                "solver_checkpoint_receipt.transport",
                format!(
                    "expected {SOLVER_CHECKPOINT_RECEIPT_TRANSPORT_BYTES} bytes, got {}",
                    bytes.len()
                ),
            ));
        }
        let u32_at = |offset: usize| {
            u32::from_le_bytes(
                bytes[offset..offset + 4]
                    .try_into()
                    .expect("transport length"),
            )
        };
        let u64_at = |offset: usize| {
            u64::from_le_bytes(
                bytes[offset..offset + 8]
                    .try_into()
                    .expect("transport length"),
            )
        };
        let hash_at = |offset: usize| {
            ContentHash(
                bytes[offset..offset + 32]
                    .try_into()
                    .expect("transport length"),
            )
        };
        let version = u32_at(0);
        if version != SOLVER_CHECKPOINT_RECEIPT_IDENTITY_VERSION {
            return Err(invalid(
                "solver_checkpoint_receipt.version",
                format!("found v{version}, require v{SOLVER_CHECKPOINT_RECEIPT_IDENTITY_VERSION}"),
            ));
        }
        let receipt = Self {
            ledger_instance_id: LedgerInstanceId(
                bytes[4..20].try_into().expect("transport length"),
            ),
            session: u64_at(20),
            run: fs_exec::RunId(u64_at(28)),
            pause_authority: hash_at(36),
            gate_generation: u64_at(68),
            solver_state_artifact: hash_at(76),
            drain_report_hash: hash_at(108),
            registered_workers: u64_at(140),
            drained_workers: u64_at(148),
            content_hash: hash_at(156),
        };
        validate_checkpoint_receipt(receipt)?;
        Ok(receipt)
    }
}

/// Canonical immutable mutation claim offered before caller work begins.
///
/// The caller supplies the opaque request authority, but fs-ledger binds it to
/// the currently checked physical ledger, governor, session-open authority,
/// exact payload bytes, scope, generation, and optional causal ordinal.
#[derive(Debug, Clone, Copy)]
pub struct SessionMutationClaim<'a> {
    /// Opaque pre-execution request authority.
    pub authority: ContentHash,
    /// Physical ledger expected by the caller.
    pub ledger_instance_id: LedgerInstanceId,
    /// Durable governor identity hash.
    pub governor_hash: ContentHash,
    /// Durable session-open authority or receipt hash.
    pub session_open_hash: ContentHash,
    /// Bounded mutation kind, for example meter or submission.
    pub kind: &'a str,
    /// Numeric session identity, stored as exact big-endian bytes.
    pub session: u64,
    /// Immutable canonical ledger namespace.
    pub ledger_scope: &'a str,
    /// Session execution generation, stored as exact big-endian bytes.
    pub generation: u64,
    /// Causal meter or mutation ordinal in `1..=i64::MAX` when present.
    /// New submission claims must carry their admission ordinal; typed reads
    /// retain compatibility with immutable v6 submission rows that stored NULL.
    pub causal_ordinal: Option<u64>,
    /// Canonical typed request payload bytes.
    pub payload: &'a [u8],
}

/// Positive authority returned only to the caller that durably inserted a
/// fresh Pending claim.
///
/// Its fields are private so an observer of Pending cannot construct permission
/// to terminalize work it did not claim.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SessionClaimPermit {
    authority: ContentHash,
    ledger_instance_id: LedgerInstanceId,
    claim_hash: ContentHash,
}

impl SessionClaimPermit {
    /// Opaque request authority owned by this permit.
    #[must_use]
    pub const fn authority(self) -> ContentHash {
        self.authority
    }
}

/// One terminal receipt offered by the owner of a fresh mutation claim.
#[derive(Debug, Clone, Copy)]
pub struct SessionTerminalRow<'a> {
    /// Exact claim whose work reached a terminal state.
    pub claim: SessionMutationClaim<'a>,
    /// Positive permit returned by a pre-execution claim insertion.
    ///
    /// None is reserved for already-completed non-execution mutations whose
    /// claim and terminal may be inserted together in the terminal batch.
    pub permit: Option<SessionClaimPermit>,
    /// Canonical typed terminal receipt bytes.
    pub receipt: &'a [u8],
}

/// One terminal receipt and the ordered global audit events it owns.
#[derive(Debug)]
pub struct SessionTerminalGroup<'a> {
    /// Immutable terminal row.
    pub terminal: SessionTerminalRow<'a>,
    /// Ordered events appended and linked atomically with the terminal.
    pub events: &'a [EventRow<'a>],
}

/// One deterministic atomic terminal batch.
///
/// The batch identity is computed internally from the checked physical ledger
/// and the complete ordered group preimage.
#[derive(Debug)]
pub struct SessionTerminalBatch<'a> {
    /// Ordered terminal/event groups in the prepared flush.
    pub groups: &'a [SessionTerminalGroup<'a>],
}

/// One bounded, hash-verified durable mutation claim.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StoredSessionMutationClaim {
    /// Opaque request authority.
    pub authority: ContentHash,
    /// Physical ledger identity recorded with the claim.
    pub ledger_instance_id: LedgerInstanceId,
    /// Durable governor identity hash.
    pub governor_hash: ContentHash,
    /// Durable session-open authority or receipt hash.
    pub session_open_hash: ContentHash,
    /// Registry row schema version.
    pub schema_version: i64,
    /// Canonical mutation kind.
    pub kind: String,
    /// Numeric session identity.
    pub session: u64,
    /// Immutable ledger scope.
    pub ledger_scope: String,
    /// Session execution generation.
    pub generation: u64,
    /// Optional causal meter or mutation ordinal.
    pub causal_ordinal: Option<u64>,
    /// Canonical typed request payload bytes.
    pub payload: Vec<u8>,
    /// Plain BLAKE3 hash recomputed over payload on every read.
    pub payload_hash: ContentHash,
    /// Domain-separated hash of the complete claim envelope.
    pub claim_hash: ContentHash,
    /// Ledger wall-clock timestamp assigned at claim insertion.
    pub created_at: i64,
}

/// One bounded, hash- and event-verified terminal receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StoredSessionTerminal {
    /// Exact durable claim terminalized by this receipt.
    pub claim: StoredSessionMutationClaim,
    /// Canonical typed terminal receipt bytes.
    pub receipt: Vec<u8>,
    /// Plain BLAKE3 hash recomputed over receipt on every read.
    pub receipt_hash: ContentHash,
    /// Number of immutable authority/sequence event links.
    pub event_count: usize,
    /// Hash of the complete ordered linked global-event group.
    pub events_hash: ContentHash,
    /// Conservatively encoded claim, terminal, and event bytes.
    pub encoded_bytes: usize,
    /// Ledger wall-clock timestamp assigned at terminal insertion.
    pub created_at: i64,
}

/// One bounded durable canonical batch marker.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StoredSessionFlushBatch {
    /// Internally computed canonical batch identity.
    pub batch_id: ContentHash,
    /// Physical ledger identity bound into the batch identity.
    pub ledger_instance_id: LedgerInstanceId,
    /// Registry row schema version.
    pub schema_version: i64,
    /// Number of ordered terminal groups.
    pub terminal_count: usize,
    /// Total owned audit events.
    pub event_count: usize,
    /// Conservatively encoded bytes in the complete batch.
    pub encoded_bytes: usize,
    /// Ledger wall-clock timestamp assigned at insertion.
    pub created_at: i64,
}

/// Result of atomically claiming one mutation before caller work begins.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SessionMutationClaimResult {
    /// This call durably created the Pending claim and alone may execute it.
    Claimed {
        /// Positive terminalization authority.
        permit: SessionClaimPermit,
    },
    /// The identical claim already exists without a terminal receipt.
    ///
    /// The caller must report Pending or Indeterminate and must not execute.
    Pending {
        /// Exact verified durable claim.
        claim: Box<StoredSessionMutationClaim>,
    },
    /// The identical claim already has a verified terminal receipt.
    Terminal {
        /// Original exact receipt and verified owned-event commitment.
        terminal: Box<StoredSessionTerminal>,
    },
}

/// Result of one atomic terminal batch write.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SessionTerminalBatchResult {
    /// The canonical batch marker was new.
    Committed {
        /// Internally computed canonical identity.
        batch_id: ContentHash,
        /// Terminal rows newly inserted by this transaction.
        terminals_inserted: usize,
        /// Global audit-event rows newly appended by this transaction.
        events_appended: usize,
    },
    /// The exact canonical batch already committed. No row was written.
    Replayed {
        /// Internally computed canonical identity.
        batch_id: ContentHash,
    },
}

#[derive(Debug, Clone, Copy)]
struct PreparedClaim {
    payload_hash: ContentHash,
    claim_hash: ContentHash,
    encoded_bytes: usize,
}

#[derive(Debug, Clone, Copy)]
struct PreparedTerminal {
    claim: PreparedClaim,
    receipt_hash: ContentHash,
    event_count: usize,
    events_hash: ContentHash,
    encoded_bytes: usize,
}

#[derive(Debug)]
struct PreparedBatch {
    ledger_instance_id: LedgerInstanceId,
    batch_id: ContentHash,
    terminals: Vec<PreparedTerminal>,
    event_count: usize,
    encoded_bytes: usize,
}

#[derive(Debug)]
struct SimpleGenerationProbe {
    claim: StoredSessionMutationClaim,
    terminalized: bool,
}

fn invalid(field: &str, problem: impl Into<String>) -> LedgerError {
    LedgerError::Invalid {
        field: field.to_string(),
        problem: problem.into(),
    }
}

fn solver_checkpoint_receipt_hash(receipt: SolverCheckpointReceipt) -> ContentHash {
    let mut hasher = Blake3::new();
    hasher.update(SOLVER_CHECKPOINT_RECEIPT_IDENTITY_DOMAIN);
    hasher.update(&SOLVER_CHECKPOINT_RECEIPT_IDENTITY_VERSION.to_le_bytes());
    hasher.update(&receipt.ledger_instance_id.as_bytes());
    hasher.update(&receipt.session.to_le_bytes());
    hasher.update(&receipt.run.0.to_le_bytes());
    hasher.update(receipt.pause_authority.as_bytes());
    hasher.update(&receipt.gate_generation.to_le_bytes());
    hasher.update(receipt.solver_state_artifact.as_bytes());
    hasher.update(receipt.drain_report_hash.as_bytes());
    hasher.update(&receipt.registered_workers.to_le_bytes());
    hasher.update(&receipt.drained_workers.to_le_bytes());
    hasher.finalize()
}

fn validate_checkpoint_receipt(receipt: SolverCheckpointReceipt) -> Result<(), LedgerError> {
    if receipt.registered_workers == 0 {
        return Err(invalid(
            "solver_checkpoint_receipt.registered_workers",
            "a zero-worker report cannot attest an executor drain",
        ));
    }
    if receipt.drained_workers != receipt.registered_workers {
        return Err(invalid(
            "solver_checkpoint_receipt.drained_workers",
            format!(
                "{} drained workers differs from {} registered workers",
                receipt.drained_workers, receipt.registered_workers
            ),
        ));
    }
    let computed = solver_checkpoint_receipt_hash(receipt);
    if computed != receipt.content_hash {
        return Err(invalid(
            "solver_checkpoint_receipt.content_hash",
            format!(
                "declared {} differs from recomputed {computed}",
                receipt.content_hash
            ),
        ));
    }
    Ok(())
}

fn stored_corrupt(authority: ContentHash, detail: impl Into<String>) -> LedgerError {
    LedgerError::Corrupt {
        hash_hex: authority.to_hex(),
        detail: format!("session mutation registry: {}", detail.into()),
    }
}

fn checked_add(current: usize, added: usize, resource: &'static str) -> Result<usize, LedgerError> {
    current.checked_add(added).ok_or_else(|| {
        invalid(
            resource,
            format!(
                "encoded byte count overflowed usize; limit is {MAX_SESSION_FLUSH_ENCODED_BYTES}"
            ),
        )
    })
}

fn require_bounded_ascii(
    field: &'static str,
    value: &str,
    max_bytes: usize,
) -> Result<(), LedgerError> {
    if value.is_empty() {
        return Err(invalid(field, "must not be empty"));
    }
    if value.len() > max_bytes {
        return Err(invalid(
            field,
            format!("{} bytes exceeds the {max_bytes}-byte limit", value.len()),
        ));
    }
    if !value.bytes().all(|byte| (b'!'..=b'~').contains(&byte)) {
        return Err(invalid(
            field,
            "must contain only canonical visible ASCII bytes",
        ));
    }
    Ok(())
}

fn update_len(hasher: &mut Blake3, len: usize) {
    hasher.update(
        &u64::try_from(len)
            .expect("bounded session-registry length fits u64")
            .to_le_bytes(),
    );
}

fn update_bytes(hasher: &mut Blake3, bytes: &[u8]) {
    update_len(hasher, bytes.len());
    hasher.update(bytes);
}

fn update_optional_u64(hasher: &mut Blake3, value: Option<u64>) {
    match value {
        Some(value) => {
            hasher.update(&[1]);
            hasher.update(&value.to_be_bytes());
        }
        None => hasher.update(&[0]),
    }
}

fn update_event_preimage(
    hasher: &mut Blake3,
    session: &[u8],
    t: i64,
    kind: &str,
    payload: Option<&str>,
) {
    update_bytes(hasher, session);
    hasher.update(&t.to_le_bytes());
    update_bytes(hasher, kind.as_bytes());
    match payload {
        Some(payload) => {
            hasher.update(&[1]);
            update_bytes(hasher, payload.as_bytes());
        }
        None => hasher.update(&[0]),
    }
}

fn events_hasher(event_count: usize) -> Blake3 {
    let mut hasher = Blake3::new();
    hasher.update(SESSION_EVENTS_HASH_DOMAIN);
    update_len(&mut hasher, event_count);
    hasher
}

fn hash_events(events: &[EventRow<'_>]) -> ContentHash {
    let mut hasher = events_hasher(events.len());
    for event in events {
        update_event_preimage(
            &mut hasher,
            event.session.unwrap_or_default(),
            event.t,
            event.kind,
            event.payload,
        );
    }
    hasher.finalize()
}

/// Compute the canonical ordered event-group hash used by durable session
/// terminals. Higher layers use this to verify that typed receipt semantics
/// reproduce the exact audit rows authenticated by [`StoredSessionTerminal`].
/// Input validation and bounds remain the responsibility of the terminal
/// batch API; this helper is a pure deterministic identity function.
#[must_use]
pub fn session_terminal_events_hash(events: &[EventRow<'_>]) -> ContentHash {
    hash_events(events)
}

fn compute_claim_hash(claim: SessionMutationClaim<'_>, payload_hash: ContentHash) -> ContentHash {
    let mut hasher = Blake3::new();
    hasher.update(SESSION_CLAIM_HASH_DOMAIN);
    hasher.update(claim.authority.as_bytes());
    hasher.update(&claim.ledger_instance_id.as_bytes());
    hasher.update(claim.governor_hash.as_bytes());
    hasher.update(claim.session_open_hash.as_bytes());
    hasher.update(&SESSION_REGISTRY_ROW_SCHEMA_VERSION.to_le_bytes());
    update_bytes(&mut hasher, claim.kind.as_bytes());
    hasher.update(&claim.session.to_be_bytes());
    update_bytes(&mut hasher, claim.ledger_scope.as_bytes());
    hasher.update(&claim.generation.to_be_bytes());
    update_optional_u64(&mut hasher, claim.causal_ordinal);
    hasher.update(payload_hash.as_bytes());
    hasher.finalize()
}

fn validate_claim(
    claim: SessionMutationClaim<'_>,
    current_ledger: LedgerInstanceId,
) -> Result<PreparedClaim, LedgerError> {
    if claim.ledger_instance_id != current_ledger {
        return Err(invalid(
            "session_claim.ledger_instance_id",
            format!(
                "claim is bound to ledger {} but this checked ledger is {}",
                claim.ledger_instance_id, current_ledger
            ),
        ));
    }
    require_bounded_ascii(
        "session_claim.kind",
        claim.kind,
        MAX_SESSION_TERMINAL_KIND_BYTES,
    )?;
    require_bounded_ascii(
        "session_claim.ledger_scope",
        claim.ledger_scope,
        MAX_SESSION_TERMINAL_SCOPE_BYTES,
    )?;
    if let Some(ordinal) = claim.causal_ordinal
        && (ordinal == 0 || ordinal > i64::MAX as u64)
    {
        return Err(invalid(
            "session_claim.causal_ordinal",
            format!("ordinal {ordinal} is outside the canonical 1..=i64::MAX range"),
        ));
    }
    if claim.kind == PRECLAIM_REQUIRED_SUBMISSION_KIND && claim.causal_ordinal.is_none() {
        return Err(invalid(
            "session_claim.causal_ordinal",
            "pre-execution submission claims must bind their admission ordinal",
        ));
    }
    if claim.payload.len() > MAX_SESSION_CLAIM_PAYLOAD_BYTES {
        return Err(invalid(
            "session_claim.payload",
            format!(
                "{} bytes exceeds the {MAX_SESSION_CLAIM_PAYLOAD_BYTES}-byte limit",
                claim.payload.len()
            ),
        ));
    }
    let encoded_bytes = SESSION_CLAIM_ROW_FRAMING_BYTES
        .checked_add(claim.kind.len())
        .and_then(|bytes| bytes.checked_add(claim.ledger_scope.len()))
        .and_then(|bytes| bytes.checked_add(claim.payload.len()))
        .ok_or_else(|| {
            invalid(
                "session_flush_encoded_bytes",
                "claim encoded byte count overflowed usize",
            )
        })?;
    let payload_hash = hash_bytes(claim.payload);
    Ok(PreparedClaim {
        payload_hash,
        claim_hash: compute_claim_hash(claim, payload_hash),
        encoded_bytes,
    })
}

fn validate_event(
    ledger: &Ledger,
    claim: SessionMutationClaim<'_>,
    event: &EventRow<'_>,
) -> Result<usize, LedgerError> {
    let session = claim.session.to_be_bytes();
    if event.session != Some(session.as_slice()) {
        return Err(invalid(
            "session_terminal_event.session",
            format!(
                "owned event must carry the exact 8-byte session {} identity",
                claim.session
            ),
        ));
    }
    require_bounded_ascii(
        "session_terminal_event.kind",
        event.kind,
        MAX_SESSION_TERMINAL_EVENT_KIND_BYTES,
    )?;
    if event.kind == VCS_IDENTITY_EVENT_KIND {
        return Err(invalid(
            "session_terminal_event.kind",
            format!("event kind {VCS_IDENTITY_EVENT_KIND:?} is reserved"),
        ));
    }
    let payload_bytes = event.payload.map_or(0, str::len);
    if payload_bytes > MAX_SESSION_TERMINAL_EVENT_PAYLOAD_BYTES {
        return Err(invalid(
            "session_terminal_event.payload",
            format!(
                "{payload_bytes} bytes exceeds the \
                 {MAX_SESSION_TERMINAL_EVENT_PAYLOAD_BYTES}-byte limit"
            ),
        ));
    }
    if let Some(payload) = event.payload {
        ledger.require_json("session_terminal_event.payload", payload, false)?;
    }
    SESSION_EVENT_ROW_FRAMING_BYTES
        .checked_add(session.len())
        .and_then(|bytes| bytes.checked_add(event.kind.len()))
        .and_then(|bytes| bytes.checked_add(payload_bytes))
        .ok_or_else(|| {
            invalid(
                "session_flush_encoded_bytes",
                "event encoded byte count overflowed usize",
            )
        })
}

fn validate_terminal(
    ledger: &Ledger,
    group: &SessionTerminalGroup<'_>,
    current_ledger: LedgerInstanceId,
) -> Result<PreparedTerminal, LedgerError> {
    let row = group.terminal;
    let claim = validate_claim(row.claim, current_ledger)?;
    if let Some(permit) = row.permit
        && (permit.authority != row.claim.authority
            || permit.ledger_instance_id != current_ledger
            || permit.claim_hash != claim.claim_hash)
    {
        return Err(invalid(
            "session_terminal.permit",
            "terminalization permit does not match the exact durable claim",
        ));
    }
    if row.receipt.is_empty() || row.receipt.len() > MAX_SESSION_TERMINAL_RECEIPT_BYTES {
        return Err(invalid(
            "session_terminal.receipt",
            format!(
                "{} bytes is outside the 1..={MAX_SESSION_TERMINAL_RECEIPT_BYTES} bound",
                row.receipt.len()
            ),
        ));
    }
    if group.events.len() > MAX_SESSION_FLUSH_EVENTS {
        return Err(invalid(
            "session_terminal.events",
            format!(
                "{} rows exceeds the {MAX_SESSION_FLUSH_EVENTS}-event batch limit",
                group.events.len()
            ),
        ));
    }

    let mut encoded_bytes = claim
        .encoded_bytes
        .checked_add(SESSION_TERMINAL_ROW_FRAMING_BYTES)
        .and_then(|bytes| bytes.checked_add(row.receipt.len()))
        .ok_or_else(|| {
            invalid(
                "session_flush_encoded_bytes",
                "terminal encoded byte count overflowed usize",
            )
        })?;
    for event in group.events {
        encoded_bytes = checked_add(
            encoded_bytes,
            validate_event(ledger, row.claim, event)?,
            "session_flush_encoded_bytes",
        )?;
    }
    if encoded_bytes > MAX_SESSION_FLUSH_ENCODED_BYTES {
        return Err(invalid(
            "session_flush_encoded_bytes",
            format!(
                "{encoded_bytes} bytes exceeds the \
                 {MAX_SESSION_FLUSH_ENCODED_BYTES}-byte limit"
            ),
        ));
    }

    Ok(PreparedTerminal {
        claim,
        receipt_hash: hash_bytes(row.receipt),
        event_count: group.events.len(),
        events_hash: hash_events(group.events),
        encoded_bytes,
    })
}

fn update_prepared_terminal_preimage(
    hasher: &mut Blake3,
    authority: ContentHash,
    prepared: PreparedTerminal,
) {
    hasher.update(authority.as_bytes());
    hasher.update(prepared.claim.claim_hash.as_bytes());
    hasher.update(prepared.receipt_hash.as_bytes());
    update_len(hasher, prepared.event_count);
    hasher.update(prepared.events_hash.as_bytes());
    update_len(hasher, prepared.encoded_bytes);
}

fn prepare_batch(
    ledger: &Ledger,
    batch: &SessionTerminalBatch<'_>,
) -> Result<PreparedBatch, LedgerError> {
    if batch.groups.is_empty() || batch.groups.len() > MAX_SESSION_FLUSH_TERMINALS {
        return Err(invalid(
            "session_flush_batch.groups",
            format!(
                "{} terminal rows is outside the 1..={MAX_SESSION_FLUSH_TERMINALS} bound",
                batch.groups.len()
            ),
        ));
    }
    let current_ledger = ledger.checked_instance_id()?;
    let mut authorities = BTreeSet::new();
    let mut terminals = Vec::with_capacity(batch.groups.len());
    let mut event_count = 0usize;
    let mut encoded_bytes = 0usize;
    for group in batch.groups {
        if !authorities.insert(group.terminal.claim.authority) {
            return Err(invalid(
                "session_flush_batch.groups",
                format!(
                    "claim authority {} appears more than once in one batch",
                    group.terminal.claim.authority
                ),
            ));
        }
        let prepared = validate_terminal(ledger, group, current_ledger)?;
        event_count = checked_add(
            event_count,
            prepared.event_count,
            "session_flush_event_count",
        )?;
        if event_count > MAX_SESSION_FLUSH_EVENTS {
            return Err(invalid(
                "session_flush_event_count",
                format!("{event_count} rows exceeds the {MAX_SESSION_FLUSH_EVENTS}-event limit"),
            ));
        }
        encoded_bytes = checked_add(
            encoded_bytes,
            prepared.encoded_bytes,
            "session_flush_encoded_bytes",
        )?;
        if encoded_bytes > MAX_SESSION_FLUSH_ENCODED_BYTES {
            return Err(invalid(
                "session_flush_encoded_bytes",
                format!(
                    "{encoded_bytes} bytes exceeds the \
                     {MAX_SESSION_FLUSH_ENCODED_BYTES}-byte limit"
                ),
            ));
        }
        terminals.push(prepared);
    }

    let mut hasher = Blake3::new();
    hasher.update(SESSION_BATCH_HASH_DOMAIN);
    hasher.update(&current_ledger.as_bytes());
    hasher.update(&SESSION_REGISTRY_ROW_SCHEMA_VERSION.to_le_bytes());
    update_len(&mut hasher, batch.groups.len());
    for (group, prepared) in batch.groups.iter().zip(&terminals) {
        update_prepared_terminal_preimage(&mut hasher, group.terminal.claim.authority, *prepared);
    }
    Ok(PreparedBatch {
        ledger_instance_id: current_ledger,
        batch_id: hasher.finalize(),
        terminals,
        event_count,
        encoded_bytes,
    })
}

fn row_blob<const N: usize>(
    row: &fsqlite::Row,
    index: usize,
    authority: ContentHash,
    field: &'static str,
) -> Result<[u8; N], LedgerError> {
    let Some(SqliteValue::Blob(bytes)) = row.get(index) else {
        return Err(stored_corrupt(authority, format!("{field} is not a BLOB")));
    };
    bytes
        .as_ref()
        .try_into()
        .map_err(|_| stored_corrupt(authority, format!("{field} is not exactly {N} bytes")))
}

fn row_blob_vec(
    row: &fsqlite::Row,
    index: usize,
    authority: ContentHash,
    field: &'static str,
) -> Result<Vec<u8>, LedgerError> {
    let Some(SqliteValue::Blob(bytes)) = row.get(index) else {
        return Err(stored_corrupt(authority, format!("{field} is not a BLOB")));
    };
    Ok(bytes.as_ref().to_vec())
}

fn row_optional_u64(
    row: &fsqlite::Row,
    index: usize,
    authority: ContentHash,
    field: &'static str,
) -> Result<Option<u64>, LedgerError> {
    match row.get(index) {
        Some(SqliteValue::Null) => Ok(None),
        Some(SqliteValue::Blob(_)) => Ok(Some(u64::from_be_bytes(row_blob(
            row, index, authority, field,
        )?))),
        _ => Err(stored_corrupt(
            authority,
            format!("{field} is neither NULL nor an 8-byte BLOB"),
        )),
    }
}

fn row_usize(
    row: &fsqlite::Row,
    index: usize,
    authority: ContentHash,
    field: &'static str,
) -> Result<usize, LedgerError> {
    let value = row_i64(row, index, "session registry bounded read")
        .map_err(|_| stored_corrupt(authority, format!("{field} is not an INTEGER")))?;
    usize::try_from(value)
        .map_err(|_| stored_corrupt(authority, format!("{field} is negative or too large")))
}

fn row_i64_registry(
    row: &fsqlite::Row,
    index: usize,
    authority: ContentHash,
    field: &'static str,
) -> Result<i64, LedgerError> {
    row_i64(row, index, "session registry bounded read")
        .map_err(|_| stored_corrupt(authority, format!("{field} is not an INTEGER")))
}

fn decode_stored_session_claim(
    row: &fsqlite::Row,
    offset: usize,
    authority: ContentHash,
    current_ledger: LedgerInstanceId,
) -> Result<StoredSessionMutationClaim, LedgerError> {
    let stored_ledger: [u8; 16] =
        row_blob(row, offset, authority, "session_claim.ledger_instance_id")?;
    if stored_ledger != current_ledger.as_bytes() {
        return Err(stored_corrupt(
            authority,
            format!(
                "claim belongs to ledger bytes {}, not current ledger {current_ledger}",
                hex_bytes(&stored_ledger)
            ),
        ));
    }

    let payload = row_blob_vec(row, offset + 9, authority, "payload")?;
    let payload_hash = ContentHash(row_blob(row, offset + 10, authority, "payload_hash")?);
    if hash_bytes(&payload) != payload_hash {
        return Err(stored_corrupt(
            authority,
            "payload bytes do not match payload_hash",
        ));
    }
    let claim_hash = ContentHash(row_blob(row, offset + 11, authority, "claim_hash")?);
    let stored = StoredSessionMutationClaim {
        authority,
        ledger_instance_id: current_ledger,
        governor_hash: ContentHash(row_blob(row, offset + 1, authority, "governor_hash")?),
        session_open_hash: ContentHash(row_blob(row, offset + 2, authority, "session_open_hash")?),
        schema_version: row_i64_registry(row, offset + 3, authority, "registry_schema_version")?,
        kind: row_text(row, offset + 4, "session claim kind")?,
        session: u64::from_be_bytes(row_blob(row, offset + 5, authority, "session")?),
        ledger_scope: row_text(row, offset + 6, "session claim ledger_scope")?,
        generation: u64::from_be_bytes(row_blob(row, offset + 7, authority, "generation")?),
        causal_ordinal: row_optional_u64(row, offset + 8, authority, "causal_ordinal")?,
        payload,
        payload_hash,
        claim_hash,
        created_at: row_i64_registry(row, offset + 12, authority, "created_at")?,
    };
    if compute_claim_hash(stored.as_input(), stored.payload_hash) != stored.claim_hash {
        return Err(stored_corrupt(
            authority,
            "claim envelope does not match claim_hash",
        ));
    }
    Ok(stored)
}

impl StoredSessionMutationClaim {
    fn as_input(&self) -> SessionMutationClaim<'_> {
        SessionMutationClaim {
            authority: self.authority,
            ledger_instance_id: self.ledger_instance_id,
            governor_hash: self.governor_hash,
            session_open_hash: self.session_open_hash,
            kind: &self.kind,
            session: self.session,
            ledger_scope: &self.ledger_scope,
            generation: self.generation,
            causal_ordinal: self.causal_ordinal,
            payload: &self.payload,
        }
    }

    fn matches(&self, offered: SessionMutationClaim<'_>, prepared: PreparedClaim) -> bool {
        self.authority == offered.authority
            && self.ledger_instance_id == offered.ledger_instance_id
            && self.governor_hash == offered.governor_hash
            && self.session_open_hash == offered.session_open_hash
            && self.schema_version == SESSION_REGISTRY_ROW_SCHEMA_VERSION
            && self.kind == offered.kind
            && self.session == offered.session
            && self.ledger_scope == offered.ledger_scope
            && self.generation == offered.generation
            && self.causal_ordinal == offered.causal_ordinal
            && self.payload.as_slice() == offered.payload
            && self.payload_hash == prepared.payload_hash
            && self.claim_hash == prepared.claim_hash
    }
}

impl StoredSessionTerminal {
    fn matches(&self, offered: SessionTerminalRow<'_>, prepared: PreparedTerminal) -> bool {
        self.claim.matches(offered.claim, prepared.claim)
            && self.receipt.as_slice() == offered.receipt
            && self.receipt_hash == prepared.receipt_hash
            && self.event_count == prepared.event_count
            && self.events_hash == prepared.events_hash
            && self.encoded_bytes == prepared.encoded_bytes
    }
}

impl StoredSessionFlushBatch {
    fn matches(&self, prepared: &PreparedBatch, terminal_count: usize) -> bool {
        self.batch_id == prepared.batch_id
            && self.ledger_instance_id == prepared.ledger_instance_id
            && self.schema_version == SESSION_REGISTRY_ROW_SCHEMA_VERSION
            && self.terminal_count == terminal_count
            && self.event_count == prepared.event_count
            && self.encoded_bytes == prepared.encoded_bytes
    }
}

impl Ledger {
    fn validate_solver_checkpoint_artifact(
        &self,
        artifact: ContentHash,
        run: fs_exec::RunId,
    ) -> Result<(), LedgerError> {
        let info = self.artifact_info(&artifact)?.ok_or_else(|| {
            invalid(
                "solver_checkpoint_receipt.solver_state_artifact",
                format!("artifact {artifact} does not exist"),
            )
        })?;
        if info.kind != SOLVER_STATE_ARTIFACT_KIND {
            return Err(invalid(
                "solver_checkpoint_receipt.solver_state_artifact",
                format!(
                    "artifact {artifact} has kind {:?}, require {SOLVER_STATE_ARTIFACT_KIND:?}",
                    info.kind
                ),
            ));
        }
        let bytes = self
            .get_artifact_bounded(&artifact, MAX_SOLVER_CHECKPOINT_ARTIFACT_BYTES)?
            .ok_or_else(|| {
                invalid(
                    "solver_checkpoint_receipt.solver_state_artifact",
                    format!("artifact {artifact} disappeared during validation"),
                )
            })?;
        let envelope = fs_exec::solver::envelope::inspect(&bytes).map_err(|error| {
            invalid(
                "solver_checkpoint_receipt.solver_state_artifact",
                format!("artifact {artifact} is not a valid solver-state envelope: {error}"),
            )
        })?;
        if envelope.provenance() != run.0 {
            return Err(invalid(
                "solver_checkpoint_receipt.run",
                format!(
                    "snapshot provenance {} differs from drained executor run {}",
                    envelope.provenance(),
                    run.0
                ),
            ));
        }
        Ok(())
    }

    fn solver_checkpoint_receipt_at_instance(
        &self,
        pause_authority: ContentHash,
        current_ledger: LedgerInstanceId,
    ) -> Result<Option<SolverCheckpointReceipt>, LedgerError> {
        let rows = self
            .conn
            .query_with_params(
                "SELECT receipt_hash, ledger_instance_id, session, run, pause_authority, \
                        gate_generation, solver_state_artifact, drain_report_hash, \
                        registered_workers, drained_workers \
                 FROM session_checkpoint_receipts \
                 WHERE pause_authority = ?1 LIMIT 2",
                &[blob_param(pause_authority.as_bytes())],
            )
            .map_err(|error| sql_err("solver checkpoint receipt read", &error))?;
        if rows.is_empty() {
            return Ok(None);
        }
        if rows.len() != 1 {
            return Err(stored_corrupt(
                pause_authority,
                "solver checkpoint pause authority names multiple receipts",
            ));
        }
        let row = &rows[0];
        let receipt = SolverCheckpointReceipt {
            content_hash: ContentHash(row_blob(
                row,
                0,
                pause_authority,
                "solver_checkpoint.receipt_hash",
            )?),
            ledger_instance_id: LedgerInstanceId(row_blob(
                row,
                1,
                pause_authority,
                "solver_checkpoint.ledger_instance_id",
            )?),
            session: u64::from_be_bytes(row_blob(
                row,
                2,
                pause_authority,
                "solver_checkpoint.session",
            )?),
            run: fs_exec::RunId(u64::from_be_bytes(row_blob(
                row,
                3,
                pause_authority,
                "solver_checkpoint.run",
            )?)),
            pause_authority: ContentHash(row_blob(
                row,
                4,
                pause_authority,
                "solver_checkpoint.pause_authority",
            )?),
            gate_generation: u64::from_be_bytes(row_blob(
                row,
                5,
                pause_authority,
                "solver_checkpoint.gate_generation",
            )?),
            solver_state_artifact: ContentHash(row_blob(
                row,
                6,
                pause_authority,
                "solver_checkpoint.solver_state_artifact",
            )?),
            drain_report_hash: ContentHash(row_blob(
                row,
                7,
                pause_authority,
                "solver_checkpoint.drain_report_hash",
            )?),
            registered_workers: u64::from_be_bytes(row_blob(
                row,
                8,
                pause_authority,
                "solver_checkpoint.registered_workers",
            )?),
            drained_workers: u64::from_be_bytes(row_blob(
                row,
                9,
                pause_authority,
                "solver_checkpoint.drained_workers",
            )?),
        };
        if receipt.ledger_instance_id != current_ledger {
            return Err(stored_corrupt(
                pause_authority,
                format!(
                    "checkpoint receipt belongs to ledger {}, not current ledger {current_ledger}",
                    receipt.ledger_instance_id
                ),
            ));
        }
        if receipt.pause_authority != pause_authority {
            return Err(stored_corrupt(
                pause_authority,
                "checkpoint receipt row disagrees with its indexed pause authority",
            ));
        }
        validate_checkpoint_receipt(receipt).map_err(|error| {
            stored_corrupt(
                pause_authority,
                format!("solver checkpoint receipt identity is corrupt: {error}"),
            )
        })?;
        Ok(Some(receipt))
    }

    fn insert_solver_checkpoint_receipt(
        &self,
        receipt: SolverCheckpointReceipt,
    ) -> Result<(), LedgerError> {
        self.conn
            .prepare(
                "INSERT INTO session_checkpoint_receipts( \
                    receipt_hash, ledger_instance_id, session, run, pause_authority, \
                    gate_generation, solver_state_artifact, drain_report_hash, \
                    registered_workers, drained_workers, created_at \
                 ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11)",
            )
            .map_err(|error| sql_err("solver checkpoint receipt insert prepare", &error))?
            .execute_with_params(&[
                blob_param(receipt.content_hash.as_bytes()),
                blob_param(&receipt.ledger_instance_id.as_bytes()),
                blob_param(&receipt.session.to_be_bytes()),
                blob_param(&receipt.run.0.to_be_bytes()),
                blob_param(receipt.pause_authority.as_bytes()),
                blob_param(&receipt.gate_generation.to_be_bytes()),
                blob_param(receipt.solver_state_artifact.as_bytes()),
                blob_param(receipt.drain_report_hash.as_bytes()),
                blob_param(&receipt.registered_workers.to_be_bytes()),
                blob_param(&receipt.drained_workers.to_be_bytes()),
                SqliteValue::Integer(now_wall_ns()),
            ])
            .map_err(|error| sql_err("solver checkpoint receipt insert", &error))?;
        Ok(())
    }

    /// Mint or exactly replay one immutable solver-checkpoint receipt.
    ///
    /// The referenced artifact must already exist with kind `solver-state`,
    /// its complete bounded bytes must form a checksum-valid fs-exec snapshot,
    /// and its provenance must equal the opaque drain report's run. The pause
    /// authority is unique: a retry after response loss returns the same
    /// receipt, while any changed field conflicts atomically.
    ///
    /// # Errors
    /// Open transactions, missing/oversized/corrupt artifacts, run mismatch,
    /// malformed drain evidence, authority conflicts, and ledger failures.
    pub fn attest_solver_checkpoint(
        &self,
        claim: SolverCheckpointClaim<'_>,
    ) -> Result<SolverCheckpointReceipt, LedgerError> {
        if self.in_transaction() {
            return Err(invalid(
                "solver_checkpoint_receipt.transaction",
                "checkpoint attestation must own its transaction; commit or roll back first",
            ));
        }
        let current_ledger = self.checked_instance_id()?;
        let report = *claim.drain_report;
        if report.registered_workers() == 0
            || report.drained_workers() != report.registered_workers()
        {
            return Err(invalid(
                "solver_checkpoint_receipt.drain_report",
                "executor report does not prove a non-empty fully drained worker set",
            ));
        }
        let mut report_preimage = Vec::with_capacity(28);
        report_preimage
            .extend_from_slice(&fs_exec::DRAIN_FINALIZE_REPORT_IDENTITY_VERSION.to_le_bytes());
        report_preimage.extend_from_slice(&report.run().0.to_le_bytes());
        report_preimage.extend_from_slice(&report.registered_workers().to_le_bytes());
        report_preimage.extend_from_slice(&report.drained_workers().to_le_bytes());
        let computed_report_hash = fs_blake3::hash_domain(
            fs_exec::DRAIN_FINALIZE_REPORT_IDENTITY_DOMAIN,
            &report_preimage,
        );
        if computed_report_hash != report.content_hash() {
            return Err(invalid(
                "solver_checkpoint_receipt.drain_report",
                "executor report identity does not match its private fields",
            ));
        }
        self.validate_solver_checkpoint_artifact(claim.solver_state_artifact, report.run())?;
        let mut receipt = SolverCheckpointReceipt {
            ledger_instance_id: current_ledger,
            session: claim.session,
            run: report.run(),
            pause_authority: claim.pause_authority,
            gate_generation: claim.gate_generation,
            solver_state_artifact: claim.solver_state_artifact,
            drain_report_hash: report.content_hash(),
            registered_workers: report.registered_workers(),
            drained_workers: report.drained_workers(),
            content_hash: ContentHash([0; 32]),
        };
        receipt.content_hash = solver_checkpoint_receipt_hash(receipt);
        validate_checkpoint_receipt(receipt)?;

        self.begin()?;
        let write = (|| match self
            .solver_checkpoint_receipt_at_instance(claim.pause_authority, current_ledger)?
        {
            Some(stored) if stored == receipt => Ok(stored),
            Some(_) => Err(invalid(
                "solver_checkpoint_receipt.pause_authority",
                format!(
                    "pause authority {} already stores a different checkpoint receipt",
                    claim.pause_authority
                ),
            )),
            None => {
                self.insert_solver_checkpoint_receipt(receipt)?;
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

    /// Recover one typed checkpoint receipt by its exact pause authority.
    /// Absence is `Ok(None)`; a present row is fully reverified before return.
    ///
    /// # Errors
    /// Malformed, foreign-ledger, hash-mismatched, or orphaned receipt state
    /// fails closed.
    pub fn solver_checkpoint_receipt(
        &self,
        pause_authority: ContentHash,
    ) -> Result<Option<SolverCheckpointReceipt>, LedgerError> {
        let current_ledger = self.checked_instance_id()?;
        self.note_read_query();
        let Some(receipt) =
            self.solver_checkpoint_receipt_at_instance(pause_authority, current_ledger)?
        else {
            return Ok(None);
        };
        self.validate_solver_checkpoint_artifact(receipt.solver_state_artifact, receipt.run)?;
        Ok(Some(receipt))
    }

    /// Independently re-earn ledger membership, immutable row identity, and
    /// solver-state artifact/run binding for a transport receipt candidate.
    ///
    /// # Errors
    /// Self-inconsistent transport, foreign/missing/conflicting ledger rows,
    /// or a missing/corrupt/mismatched artifact is refused.
    pub fn verify_solver_checkpoint_receipt(
        &self,
        receipt: &SolverCheckpointReceipt,
    ) -> Result<(), LedgerError> {
        validate_checkpoint_receipt(*receipt)?;
        let current_ledger = self.checked_instance_id()?;
        if receipt.ledger_instance_id != current_ledger {
            return Err(invalid(
                "solver_checkpoint_receipt.ledger_instance_id",
                format!(
                    "receipt belongs to ledger {}, not current ledger {current_ledger}",
                    receipt.ledger_instance_id
                ),
            ));
        }
        self.note_read_query();
        let stored = self
            .solver_checkpoint_receipt_at_instance(receipt.pause_authority, current_ledger)?
            .ok_or_else(|| {
                invalid(
                    "solver_checkpoint_receipt.content_hash",
                    format!(
                        "receipt {} is not stored by this ledger",
                        receipt.content_hash
                    ),
                )
            })?;
        if stored != *receipt {
            return Err(invalid(
                "solver_checkpoint_receipt.content_hash",
                "stored checkpoint receipt differs from the supplied candidate",
            ));
        }
        self.validate_solver_checkpoint_artifact(receipt.solver_state_artifact, receipt.run)
    }

    /// Fetch and hash-verify one immutable mutation claim.
    ///
    /// Absence is Ok(None). A claim may still be Pending; use
    /// session_terminal or claim_session_mutation to distinguish state.
    ///
    /// # Errors
    /// Returns a fail-closed corruption error for malformed, foreign-ledger,
    /// future-envelope, or hash-mismatched rows.
    pub fn session_mutation_claim(
        &self,
        authority: &ContentHash,
    ) -> Result<Option<StoredSessionMutationClaim>, LedgerError> {
        let current_ledger = self.checked_instance_id()?;
        self.note_read_query();
        self.session_mutation_claim_at_instance(authority, current_ledger)
    }

    fn session_mutation_claim_at_instance(
        &self,
        authority: &ContentHash,
        current_ledger: LedgerInstanceId,
    ) -> Result<Option<StoredSessionMutationClaim>, LedgerError> {
        let presence = self
            .conn
            .query_with_params(
                "SELECT \
                    (SELECT COUNT(*) FROM session_claims WHERE authority = ?1), \
                    (SELECT COUNT(*) FROM session_claim_discovery WHERE authority = ?1)",
                &[blob_param(authority.as_bytes())],
            )
            .map_err(|error| sql_err("session claim/discovery presence", &error))?;
        let presence = presence.first().ok_or_else(|| {
            stored_corrupt(
                *authority,
                "claim/discovery presence query returned no aggregate row",
            )
        })?;
        let claim_count = row_i64_registry(presence, 0, *authority, "claim authority count")?;
        let discovery_count =
            row_i64_registry(presence, 1, *authority, "claim discovery authority count")?;
        if claim_count == 0 && discovery_count == 0 {
            return Ok(None);
        }
        if claim_count != 1 || discovery_count != 1 {
            return Err(stored_corrupt(
                *authority,
                format!(
                    "claim/discovery authority multiplicity differs: claim={claim_count}, discovery={discovery_count}"
                ),
            ));
        }
        let guarded_sql = format!(
            "SELECT ledger_instance_id, governor_hash, session_open_hash, \
                    registry_schema_version, kind, session, ledger_scope, generation, \
                    causal_ordinal, payload, payload_hash, claim_hash, created_at \
             FROM session_claims WHERE authority = ?1 AND \
               typeof(ledger_instance_id) = 'blob' AND length(ledger_instance_id) = 16 AND \
               typeof(governor_hash) = 'blob' AND length(governor_hash) = 32 AND \
               typeof(session_open_hash) = 'blob' AND length(session_open_hash) = 32 AND \
               typeof(registry_schema_version) = 'integer' AND \
               registry_schema_version = {SESSION_REGISTRY_ROW_SCHEMA_VERSION} AND \
               typeof(kind) = 'text' AND \
               length(CAST(kind AS BLOB)) BETWEEN 1 AND {MAX_SESSION_TERMINAL_KIND_BYTES} AND \
               length(CAST(kind AS BLOB)) = length(kind) AND kind NOT GLOB '*[^!-~]*' AND \
               typeof(session) = 'blob' AND length(session) = 8 AND \
               typeof(ledger_scope) = 'text' AND \
               length(CAST(ledger_scope AS BLOB)) BETWEEN 1 AND {MAX_SESSION_TERMINAL_SCOPE_BYTES} AND \
               length(CAST(ledger_scope AS BLOB)) = length(ledger_scope) AND \
               ledger_scope NOT GLOB '*[^!-~]*' AND \
               typeof(generation) = 'blob' AND length(generation) = 8 AND \
               (causal_ordinal IS NULL OR \
                (typeof(causal_ordinal) = 'blob' AND length(causal_ordinal) = 8 AND \
                 causal_ordinal > X'0000000000000000' AND \
                 causal_ordinal <= X'7FFFFFFFFFFFFFFF')) AND \
               typeof(payload) = 'blob' AND \
               length(payload) BETWEEN 0 AND {MAX_SESSION_CLAIM_PAYLOAD_BYTES} AND \
               typeof(payload_hash) = 'blob' AND length(payload_hash) = 32 AND \
               typeof(claim_hash) = 'blob' AND length(claim_hash) = 32 AND \
               typeof(created_at) = 'integer' LIMIT 1"
        );
        let rows = self
            .conn
            .query_with_params(&guarded_sql, &[blob_param(authority.as_bytes())])
            .map_err(|error| sql_err("session claim guarded fetch", &error))?;
        let row = rows.first().ok_or_else(|| {
            stored_corrupt(
                *authority,
                "claim violates a type, schema-version, canonical-text, or byte bound",
            )
        })?;
        let claim = decode_stored_session_claim(row, 0, *authority, current_ledger)?;
        self.verify_session_claim_discovery(&claim)?;
        Ok(Some(claim))
    }

    pub(crate) fn verify_session_claim_discovery_backfill(&self) -> Result<(), LedgerError> {
        let current_ledger = self.read_current_instance_id()?;
        let mut cursor: Option<ContentHash> = None;
        loop {
            let (sql, params) = if let Some(after) = cursor {
                (
                    format!(
                        "SELECT authority FROM session_claims WHERE authority > ?1 \
                         ORDER BY authority LIMIT {}",
                        MAX_SESSION_FLUSH_TERMINALS + 1
                    ),
                    vec![blob_param(after.as_bytes())],
                )
            } else {
                (
                    format!(
                        "SELECT authority FROM session_claims ORDER BY authority LIMIT {}",
                        MAX_SESSION_FLUSH_TERMINALS + 1
                    ),
                    Vec::new(),
                )
            };
            let rows = self
                .conn
                .query_with_params(&sql, &params)
                .map_err(|error| sql_err("session claim discovery backfill scan", &error))?;
            let has_more = rows.len() > MAX_SESSION_FLUSH_TERMINALS;
            let mut last = None;
            for row in rows.iter().take(MAX_SESSION_FLUSH_TERMINALS) {
                let authority = ContentHash(row_blob(
                    row,
                    0,
                    ContentHash([0; 32]),
                    "session claim discovery backfill authority",
                )?);
                last = Some(authority);
                self.session_mutation_claim_at_instance(&authority, current_ledger)?
                    .ok_or_else(|| {
                        stored_corrupt(
                            authority,
                            "claim disappeared during discovery backfill verification",
                        )
                    })?;
            }
            if !has_more {
                break;
            }
            cursor = Some(last.ok_or_else(|| {
                stored_corrupt(
                    ContentHash([0; 32]),
                    "claim discovery backfill page advertised a successor without a cursor",
                )
            })?);
        }
        let counts = self
            .conn
            .query(
                "SELECT \
                    (SELECT COUNT(*) FROM session_claims), \
                    (SELECT COUNT(*) FROM session_claim_discovery)",
            )
            .map_err(|error| sql_err("session claim discovery backfill counts", &error))?;
        let row = counts.first().ok_or_else(|| {
            stored_corrupt(
                ContentHash([0; 32]),
                "claim discovery backfill count query returned no aggregate row",
            )
        })?;
        let claim_count = row_i64_registry(row, 0, ContentHash([0; 32]), "backfill claim count")?;
        let discovery_count =
            row_i64_registry(row, 1, ContentHash([0; 32]), "backfill discovery count")?;
        if claim_count != discovery_count {
            return Err(stored_corrupt(
                ContentHash([0; 32]),
                format!(
                    "claim discovery backfill counts differ: claim={claim_count}, discovery={discovery_count}"
                ),
            ));
        }
        Ok(())
    }

    fn verify_session_claim_discovery(
        &self,
        claim: &StoredSessionMutationClaim,
    ) -> Result<(), LedgerError> {
        let guarded_sql = format!(
            "SELECT ledger_instance_id, governor_hash, session_open_hash, \
                    registry_schema_version, kind, session, ledger_scope, generation, \
                    causal_ordinal, payload_hash, claim_hash \
             FROM session_claim_discovery WHERE authority = ?1 AND \
               typeof(ledger_instance_id) = 'blob' AND length(ledger_instance_id) = 16 AND \
               typeof(governor_hash) = 'blob' AND length(governor_hash) = 32 AND \
               typeof(session_open_hash) = 'blob' AND length(session_open_hash) = 32 AND \
               typeof(registry_schema_version) = 'integer' AND \
               registry_schema_version = {SESSION_REGISTRY_ROW_SCHEMA_VERSION} AND \
               typeof(kind) = 'text' AND \
               length(CAST(kind AS BLOB)) BETWEEN 1 AND {MAX_SESSION_TERMINAL_KIND_BYTES} AND \
               length(CAST(kind AS BLOB)) = length(kind) AND kind NOT GLOB '*[^!-~]*' AND \
               typeof(session) = 'blob' AND length(session) = 8 AND \
               typeof(ledger_scope) = 'text' AND \
               length(CAST(ledger_scope AS BLOB)) BETWEEN 1 AND {MAX_SESSION_TERMINAL_SCOPE_BYTES} AND \
               length(CAST(ledger_scope AS BLOB)) = length(ledger_scope) AND \
               ledger_scope NOT GLOB '*[^!-~]*' AND \
               typeof(generation) = 'blob' AND length(generation) = 8 AND \
               (causal_ordinal IS NULL OR \
                (typeof(causal_ordinal) = 'blob' AND length(causal_ordinal) = 8 AND \
                 causal_ordinal > X'0000000000000000' AND \
                 causal_ordinal <= X'7FFFFFFFFFFFFFFF')) AND \
               typeof(payload_hash) = 'blob' AND length(payload_hash) = 32 AND \
               typeof(claim_hash) = 'blob' AND length(claim_hash) = 32 LIMIT 1"
        );
        let rows = self
            .conn
            .query_with_params(&guarded_sql, &[blob_param(claim.authority.as_bytes())])
            .map_err(|error| sql_err("session claim discovery guarded fetch", &error))?;
        let row = rows.first().ok_or_else(|| {
            stored_corrupt(
                claim.authority,
                "claim discovery violates a type, schema-version, canonical-text, or byte bound",
            )
        })?;
        let ledger_instance_id: [u8; 16] = row_blob(
            row,
            0,
            claim.authority,
            "claim discovery ledger_instance_id",
        )?;
        let governor_hash = ContentHash(row_blob(
            row,
            1,
            claim.authority,
            "claim discovery governor_hash",
        )?);
        let session_open_hash = ContentHash(row_blob(
            row,
            2,
            claim.authority,
            "claim discovery session_open_hash",
        )?);
        let schema_version = row_i64_registry(
            row,
            3,
            claim.authority,
            "claim discovery registry_schema_version",
        )?;
        let kind = row_text(row, 4, "claim discovery kind")?;
        let session = u64::from_be_bytes(row_blob(
            row,
            5,
            claim.authority,
            "claim discovery session",
        )?);
        let ledger_scope = row_text(row, 6, "claim discovery ledger_scope")?;
        let generation = u64::from_be_bytes(row_blob(
            row,
            7,
            claim.authority,
            "claim discovery generation",
        )?);
        let causal_ordinal =
            row_optional_u64(row, 8, claim.authority, "claim discovery causal_ordinal")?;
        let payload_hash = ContentHash(row_blob(
            row,
            9,
            claim.authority,
            "claim discovery payload_hash",
        )?);
        let claim_hash = ContentHash(row_blob(
            row,
            10,
            claim.authority,
            "claim discovery claim_hash",
        )?);
        if ledger_instance_id != claim.ledger_instance_id.as_bytes()
            || governor_hash != claim.governor_hash
            || session_open_hash != claim.session_open_hash
            || schema_version != claim.schema_version
            || kind != claim.kind
            || session != claim.session
            || ledger_scope != claim.ledger_scope
            || generation != claim.generation
            || causal_ordinal != claim.causal_ordinal
            || payload_hash != claim.payload_hash
            || claim_hash != claim.claim_hash
        {
            return Err(stored_corrupt(
                claim.authority,
                "claim discovery witness differs from the authenticated claim envelope",
            ));
        }
        Ok(())
    }

    #[allow(clippy::too_many_lines)] // One guarded SQL shape keeps every projection bounded before materialization.
    fn prepare_simple_generation_probe(&self) -> Result<PreparedStatement<'_>, LedgerError> {
        let sql = format!(
            "SELECT identity.instance_id, \
                    claim.ledger_instance_id, claim.governor_hash, \
                    claim.session_open_hash, claim.registry_schema_version, \
                    claim.kind, claim.session, claim.ledger_scope, claim.generation, \
                    claim.causal_ordinal, claim.payload, claim.payload_hash, \
                    claim.claim_hash, claim.created_at, \
                    CASE WHEN terminal.authority IS NULL THEN 0 \
                         WHEN typeof(terminal.receipt) = 'blob' AND \
                              length(terminal.receipt) BETWEEN 1 AND \
                                  {MAX_SESSION_TERMINAL_RECEIPT_BYTES} AND \
                              typeof(terminal.receipt_hash) = 'blob' AND \
                              length(terminal.receipt_hash) = 32 AND \
                              typeof(terminal.event_count) = 'integer' AND \
                              terminal.event_count BETWEEN 0 AND {MAX_SESSION_FLUSH_EVENTS} AND \
                              typeof(terminal.events_hash) = 'blob' AND \
                              length(terminal.events_hash) = 32 AND \
                              typeof(terminal.encoded_bytes) = 'integer' AND \
                              terminal.encoded_bytes BETWEEN 1 AND \
                                  {MAX_SESSION_FLUSH_ENCODED_BYTES} AND \
                              typeof(terminal.created_at) = 'integer' THEN 1 \
                         ELSE -1 END, \
                    CASE WHEN typeof(terminal.receipt) = 'blob' AND \
                                   length(terminal.receipt) BETWEEN 1 AND \
                                       {MAX_SESSION_TERMINAL_RECEIPT_BYTES} \
                         THEN terminal.receipt ELSE NULL END, \
                    CASE WHEN typeof(terminal.receipt_hash) = 'blob' AND \
                                   length(terminal.receipt_hash) = 32 \
                         THEN terminal.receipt_hash ELSE NULL END, \
                    CASE WHEN typeof(terminal.event_count) = 'integer' AND \
                                   terminal.event_count BETWEEN 0 AND \
                                       {MAX_SESSION_FLUSH_EVENTS} \
                         THEN terminal.event_count ELSE NULL END, \
                    CASE WHEN typeof(terminal.events_hash) = 'blob' AND \
                                   length(terminal.events_hash) = 32 \
                         THEN terminal.events_hash ELSE NULL END, \
                    CASE WHEN typeof(terminal.encoded_bytes) = 'integer' AND \
                                   terminal.encoded_bytes BETWEEN 1 AND \
                                       {MAX_SESSION_FLUSH_ENCODED_BYTES} \
                         THEN terminal.encoded_bytes ELSE NULL END, \
                    CASE WHEN typeof(terminal.created_at) = 'integer' \
                         THEN terminal.created_at ELSE NULL END, \
                    CASE WHEN member.authority IS NULL THEN 0 ELSE 1 END, \
                    CASE WHEN typeof(member.batch_id) = 'blob' AND \
                                   length(member.batch_id) = 32 \
                         THEN member.batch_id ELSE NULL END, \
                    CASE WHEN typeof(member.seq) = 'integer' AND \
                                   member.seq >= 0 AND \
                                   member.seq < {MAX_SESSION_FLUSH_TERMINALS} \
                         THEN member.seq ELSE NULL END, \
                    CASE WHEN batch.batch_id IS NULL THEN 0 \
                         WHEN typeof(batch.ledger_instance_id) = 'blob' AND \
                              length(batch.ledger_instance_id) = 16 AND \
                              typeof(batch.registry_schema_version) = 'integer' AND \
                              batch.registry_schema_version = \
                                  {SESSION_REGISTRY_ROW_SCHEMA_VERSION} AND \
                              typeof(batch.terminal_count) = 'integer' AND \
                              batch.terminal_count BETWEEN 1 AND \
                                  {MAX_SESSION_FLUSH_TERMINALS} AND \
                              typeof(batch.event_count) = 'integer' AND \
                              batch.event_count BETWEEN 0 AND {MAX_SESSION_FLUSH_EVENTS} AND \
                              typeof(batch.encoded_bytes) = 'integer' AND \
                              batch.encoded_bytes BETWEEN 1 AND \
                                  {MAX_SESSION_FLUSH_ENCODED_BYTES} AND \
                              typeof(batch.created_at) = 'integer' THEN 1 \
                         ELSE -1 END, \
                    CASE WHEN typeof(batch.ledger_instance_id) = 'blob' AND \
                                   length(batch.ledger_instance_id) = 16 \
                         THEN batch.ledger_instance_id ELSE NULL END, \
                    CASE WHEN typeof(batch.registry_schema_version) = 'integer' AND \
                                   batch.registry_schema_version = \
                                       {SESSION_REGISTRY_ROW_SCHEMA_VERSION} \
                         THEN batch.registry_schema_version ELSE NULL END, \
                    CASE WHEN typeof(batch.terminal_count) = 'integer' AND \
                                   batch.terminal_count BETWEEN 1 AND \
                                       {MAX_SESSION_FLUSH_TERMINALS} \
                         THEN batch.terminal_count ELSE NULL END, \
                    CASE WHEN typeof(batch.event_count) = 'integer' AND \
                                   batch.event_count BETWEEN 0 AND \
                                       {MAX_SESSION_FLUSH_EVENTS} \
                         THEN batch.event_count ELSE NULL END, \
                    CASE WHEN typeof(batch.encoded_bytes) = 'integer' AND \
                                   batch.encoded_bytes BETWEEN 1 AND \
                                       {MAX_SESSION_FLUSH_ENCODED_BYTES} \
                         THEN batch.encoded_bytes ELSE NULL END, \
                    CASE WHEN typeof(batch.created_at) = 'integer' \
                         THEN batch.created_at ELSE NULL END, \
                    CASE WHEN event_link.authority IS NULL THEN 0 ELSE 1 END, \
                    CASE WHEN member.batch_id IS NULL THEN 0 ELSE \
                         EXISTS(SELECT 1 FROM session_flush_batch_members AS other \
                                WHERE other.batch_id = member.batch_id AND \
                                      (other.seq IS NOT member.seq OR \
                                       other.authority IS NOT claim.authority) LIMIT 1) END, \
                    CASE WHEN event_link.authority IS NULL THEN 0 \
                         WHEN typeof(event_link.seq) = 'integer' AND \
                              event_link.seq >= 0 AND \
                              event_link.seq < {MAX_SESSION_FLUSH_EVENTS} AND \
                              typeof(event_link.event_id) = 'integer' AND \
                              event_link.event_id > 0 AND \
                              owned_event.id IS NOT NULL AND \
                              typeof(owned_event.id) = 'integer' AND owned_event.id > 0 AND \
                              typeof(owned_event.session) = 'blob' AND \
                              length(owned_event.session) = 8 AND \
                              typeof(owned_event.t) = 'integer' AND \
                              typeof(owned_event.kind) = 'text' AND \
                              length(CAST(owned_event.kind AS BLOB)) BETWEEN 1 AND \
                                  {MAX_SESSION_TERMINAL_EVENT_KIND_BYTES} AND \
                              length(CAST(owned_event.kind AS BLOB)) = \
                                  length(owned_event.kind) AND \
                              owned_event.kind NOT GLOB '*[^!-~]*' AND \
                              CASE WHEN owned_event.payload IS NULL THEN 1 \
                                   WHEN typeof(owned_event.payload) = 'text' AND \
                                        length(CAST(owned_event.payload AS BLOB)) BETWEEN 1 AND \
                                            {MAX_SESSION_TERMINAL_EVENT_PAYLOAD_BYTES} \
                                   THEN json_valid(owned_event.payload) ELSE 0 END = 1 \
                         THEN 1 ELSE -1 END, \
                    CASE WHEN typeof(event_link.seq) = 'integer' AND \
                                   event_link.seq >= 0 AND \
                                   event_link.seq < {MAX_SESSION_FLUSH_EVENTS} \
                         THEN event_link.seq ELSE NULL END, \
                    CASE WHEN typeof(event_link.event_id) = 'integer' AND \
                                   event_link.event_id > 0 \
                         THEN event_link.event_id ELSE NULL END, \
                    CASE WHEN typeof(owned_event.session) = 'blob' AND \
                                   length(owned_event.session) = 8 \
                         THEN owned_event.session ELSE NULL END, \
                    CASE WHEN typeof(owned_event.t) = 'integer' \
                         THEN owned_event.t ELSE NULL END, \
                    CASE WHEN typeof(owned_event.kind) = 'text' AND \
                                   length(CAST(owned_event.kind AS BLOB)) BETWEEN 1 AND \
                                       {MAX_SESSION_TERMINAL_EVENT_KIND_BYTES} AND \
                                   length(CAST(owned_event.kind AS BLOB)) = \
                                       length(owned_event.kind) AND \
                                   owned_event.kind NOT GLOB '*[^!-~]*' \
                         THEN owned_event.kind ELSE NULL END, \
                    CASE WHEN owned_event.payload IS NULL THEN NULL \
                         WHEN typeof(owned_event.payload) = 'text' AND \
                              length(CAST(owned_event.payload AS BLOB)) BETWEEN 1 AND \
                                  {MAX_SESSION_TERMINAL_EVENT_PAYLOAD_BYTES} AND \
                              json_valid(owned_event.payload) \
                         THEN owned_event.payload ELSE NULL END \
             FROM session_claims AS claim \
             JOIN session_claim_discovery AS discovery \
                  ON discovery.authority = claim.authority AND \
                     typeof(discovery.ledger_instance_id) = 'blob' AND \
                     discovery.ledger_instance_id IS claim.ledger_instance_id AND \
                     typeof(discovery.governor_hash) = 'blob' AND \
                     discovery.governor_hash IS claim.governor_hash AND \
                     typeof(discovery.session_open_hash) = 'blob' AND \
                     discovery.session_open_hash IS claim.session_open_hash AND \
                     typeof(discovery.registry_schema_version) = 'integer' AND \
                     discovery.registry_schema_version IS claim.registry_schema_version AND \
                     typeof(discovery.kind) = 'text' AND \
                     discovery.kind IS claim.kind AND \
                     typeof(discovery.session) = 'blob' AND \
                     discovery.session IS claim.session AND \
                     typeof(discovery.ledger_scope) = 'text' AND \
                     discovery.ledger_scope IS claim.ledger_scope AND \
                     typeof(discovery.generation) = 'blob' AND \
                     discovery.generation IS claim.generation AND \
                     (discovery.causal_ordinal IS NULL OR \
                      typeof(discovery.causal_ordinal) = 'blob') AND \
                     discovery.causal_ordinal IS claim.causal_ordinal AND \
                     typeof(discovery.payload_hash) = 'blob' AND \
                     discovery.payload_hash IS claim.payload_hash AND \
                     typeof(discovery.claim_hash) = 'blob' AND \
                     discovery.claim_hash IS claim.claim_hash \
             JOIN ledger_identity AS identity ON identity.singleton = 1 AND \
                  typeof(identity.singleton) = 'integer' AND \
                  typeof(identity.instance_id) = 'blob' AND \
                  length(identity.instance_id) = 16 \
             LEFT JOIN session_terminals AS terminal \
                    ON terminal.authority = claim.authority \
             LEFT JOIN session_terminal_events AS event_link \
                    ON event_link.authority = claim.authority \
             LEFT JOIN events AS owned_event ON owned_event.id = event_link.event_id \
             LEFT JOIN session_flush_batch_members AS member \
                    ON member.authority = claim.authority \
             LEFT JOIN session_flush_batches AS batch \
                    ON batch.batch_id = member.batch_id \
             WHERE claim.authority = ?1 AND \
               NOT EXISTS(SELECT 1 FROM ledger_identity AS other_identity \
                          WHERE other_identity.singleton IS NOT identity.singleton OR \
                                other_identity.instance_id IS NOT identity.instance_id \
                          LIMIT 1) AND \
               typeof(claim.ledger_instance_id) = 'blob' AND \
               length(claim.ledger_instance_id) = 16 AND \
               typeof(claim.governor_hash) = 'blob' AND \
               length(claim.governor_hash) = 32 AND \
               typeof(claim.session_open_hash) = 'blob' AND \
               length(claim.session_open_hash) = 32 AND \
               typeof(claim.registry_schema_version) = 'integer' AND \
               claim.registry_schema_version = {SESSION_REGISTRY_ROW_SCHEMA_VERSION} AND \
               typeof(claim.kind) = 'text' AND \
               length(CAST(claim.kind AS BLOB)) BETWEEN 1 AND \
                   {MAX_SESSION_TERMINAL_KIND_BYTES} AND \
               length(CAST(claim.kind AS BLOB)) = length(claim.kind) AND \
               claim.kind NOT GLOB '*[^!-~]*' AND \
               typeof(claim.session) = 'blob' AND length(claim.session) = 8 AND \
               typeof(claim.ledger_scope) = 'text' AND \
               length(CAST(claim.ledger_scope AS BLOB)) BETWEEN 1 AND \
                   {MAX_SESSION_TERMINAL_SCOPE_BYTES} AND \
               length(CAST(claim.ledger_scope AS BLOB)) = length(claim.ledger_scope) AND \
               claim.ledger_scope NOT GLOB '*[^!-~]*' AND \
               typeof(claim.generation) = 'blob' AND length(claim.generation) = 8 AND \
               (claim.causal_ordinal IS NULL OR \
                (typeof(claim.causal_ordinal) = 'blob' AND \
                 length(claim.causal_ordinal) = 8 AND \
                 claim.causal_ordinal > X'0000000000000000' AND \
                 claim.causal_ordinal <= X'7FFFFFFFFFFFFFFF')) AND \
               typeof(claim.payload) = 'blob' AND \
               length(claim.payload) BETWEEN 0 AND {MAX_SESSION_CLAIM_PAYLOAD_BYTES} AND \
               typeof(claim.payload_hash) = 'blob' AND \
               length(claim.payload_hash) = 32 AND \
               typeof(claim.claim_hash) = 'blob' AND length(claim.claim_hash) = 32 AND \
               typeof(claim.created_at) = 'integer' \
             ORDER BY member.batch_id LIMIT 2"
        );
        self.conn
            .prepare(&sql)
            .map_err(|error| sql_err("session simple generation probe prepare", &error))
    }

    /// Verify the common zero- or one-event, singleton-batch claim state with
    /// one bounded query. Any other shape returns `None` so the complete public
    /// readers can retain authority over complex and suspicious storage.
    #[allow(clippy::too_many_lines)] // Exact Rust reauthentication mirrors every bounded SQL projection.
    fn simple_generation_probe_with_statement(
        &self,
        statement: &PreparedStatement<'_>,
        authority: ContentHash,
    ) -> Result<Option<SimpleGenerationProbe>, LedgerError> {
        let rows = statement
            .query_with_params(&[blob_param(authority.as_bytes())])
            .map_err(|error| sql_err("session simple generation probe", &error))?;
        if rows.len() != 1 {
            return Ok(None);
        }
        let row = &rows[0];
        let persisted_ledger: [u8; 16] =
            row_blob(row, 0, authority, "ledger_identity.instance_id")?;
        if persisted_ledger != self.instance_id.as_bytes() {
            return Ok(None);
        }
        let claim = decode_stored_session_claim(row, 1, authority, self.instance_id)?;
        let terminal_state = row_i64_registry(row, 14, authority, "terminal state")?;
        let member_present = row_i64_registry(row, 21, authority, "batch member presence")?;
        let has_event_links = row_i64_registry(row, 31, authority, "event-link presence")?;
        if terminal_state == 0 {
            if member_present == 0 && has_event_links == 0 {
                return Ok(Some(SimpleGenerationProbe {
                    claim,
                    terminalized: false,
                }));
            }
            return Ok(None);
        }
        if terminal_state != 1 || member_present != 1 {
            return Ok(None);
        }
        if row_i64_registry(row, 32, authority, "other batch-member presence")? != 0
            || row_i64_registry(row, 24, authority, "flush-batch state")? != 1
        {
            return Ok(None);
        }

        let receipt = row_blob_vec(row, 15, authority, "receipt")?;
        let receipt_hash = ContentHash(row_blob(row, 16, authority, "receipt_hash")?);
        if hash_bytes(&receipt) != receipt_hash {
            return Ok(None);
        }
        let event_count = row_usize(row, 17, authority, "event_count")?;
        if event_count > 1 {
            return Ok(None);
        }
        let stored_events_hash = ContentHash(row_blob(row, 18, authority, "events_hash")?);
        let event_state = row_i64_registry(row, 33, authority, "owned-event state")?;
        let (events_hash, event_encoded_bytes) = match event_count {
            0 if has_event_links == 0 && event_state == 0 => (events_hasher(0).finalize(), 0),
            1 if has_event_links == 1 && event_state == 1 => {
                if row_usize(row, 34, authority, "event_link.seq")? != 0
                    || row_i64_registry(row, 35, authority, "event_link.event_id")? <= 0
                {
                    return Ok(None);
                }
                let event_session: [u8; 8] = row_blob(row, 36, authority, "event.session")?;
                if event_session != claim.session.to_be_bytes() {
                    return Ok(None);
                }
                let event_t = row_i64_registry(row, 37, authority, "event.t")?;
                let event_kind = row_text(row, 38, "session terminal event kind")?;
                let event_payload = match row.get(39) {
                    Some(SqliteValue::Null) => None,
                    Some(SqliteValue::Text(payload)) => Some(payload.as_str()),
                    _ => return Ok(None),
                };
                let mut hasher = events_hasher(1);
                update_event_preimage(
                    &mut hasher,
                    &event_session,
                    event_t,
                    &event_kind,
                    event_payload,
                );
                let event_bytes = SESSION_EVENT_ROW_FRAMING_BYTES
                    .checked_add(event_session.len())
                    .and_then(|bytes| bytes.checked_add(event_kind.len()))
                    .and_then(|bytes| bytes.checked_add(event_payload.map_or(0, str::len)))
                    .ok_or_else(|| {
                        stored_corrupt(authority, "event encoded-byte count overflowed usize")
                    })?;
                (hasher.finalize(), event_bytes)
            }
            _ => return Ok(None),
        };
        if events_hash != stored_events_hash {
            return Ok(None);
        }
        let encoded_bytes = row_usize(row, 19, authority, "encoded_bytes")?;
        let claim_encoded_bytes = SESSION_CLAIM_ROW_FRAMING_BYTES
            .checked_add(claim.kind.len())
            .and_then(|bytes| bytes.checked_add(claim.ledger_scope.len()))
            .and_then(|bytes| bytes.checked_add(claim.payload.len()))
            .ok_or_else(|| {
                stored_corrupt(authority, "claim encoded-byte count overflowed usize")
            })?;
        let expected_encoded_bytes = claim_encoded_bytes
            .checked_add(SESSION_TERMINAL_ROW_FRAMING_BYTES)
            .and_then(|bytes| bytes.checked_add(receipt.len()))
            .and_then(|bytes| bytes.checked_add(event_encoded_bytes))
            .ok_or_else(|| {
                stored_corrupt(authority, "terminal encoded-byte count overflowed usize")
            })?;
        if encoded_bytes != expected_encoded_bytes
            || encoded_bytes > MAX_SESSION_FLUSH_ENCODED_BYTES
        {
            return Ok(None);
        }
        let _ = row_i64_registry(row, 20, authority, "terminal.created_at")?;

        let batch_id = ContentHash(row_blob(row, 22, authority, "batch_member.batch_id")?);
        if row_usize(row, 23, authority, "batch_member.seq")? != 0 {
            return Ok(None);
        }
        let batch_ledger: [u8; 16] =
            row_blob(row, 25, authority, "flush_batch.ledger_instance_id")?;
        if batch_ledger != self.instance_id.as_bytes()
            || row_i64_registry(row, 26, authority, "flush_batch.registry_schema_version")?
                != SESSION_REGISTRY_ROW_SCHEMA_VERSION
            || row_usize(row, 27, authority, "flush_batch.terminal_count")? != 1
            || row_usize(row, 28, authority, "flush_batch.event_count")? != event_count
            || row_usize(row, 29, authority, "flush_batch.encoded_bytes")? != encoded_bytes
        {
            return Ok(None);
        }
        let _ = row_i64_registry(row, 30, authority, "flush_batch.created_at")?;

        let prepared = PreparedTerminal {
            claim: PreparedClaim {
                payload_hash: claim.payload_hash,
                claim_hash: claim.claim_hash,
                encoded_bytes: claim_encoded_bytes,
            },
            receipt_hash,
            event_count,
            events_hash,
            encoded_bytes,
        };
        let mut hasher = Blake3::new();
        hasher.update(SESSION_BATCH_HASH_DOMAIN);
        hasher.update(&self.instance_id.as_bytes());
        hasher.update(&SESSION_REGISTRY_ROW_SCHEMA_VERSION.to_le_bytes());
        update_len(&mut hasher, 1);
        update_prepared_terminal_preimage(&mut hasher, authority, prepared);
        if hasher.finalize() != batch_id {
            return Ok(None);
        }
        Ok(Some(SimpleGenerationProbe {
            claim,
            terminalized: true,
        }))
    }

    #[cfg(test)]
    fn simple_generation_probe(
        &self,
        authority: ContentHash,
    ) -> Result<Option<SimpleGenerationProbe>, LedgerError> {
        let statement = self.prepare_simple_generation_probe()?;
        self.simple_generation_probe_with_statement(&statement, authority)
    }

    /// Count every immutable mutation claim in one durable governor namespace.
    ///
    /// This constant-space snapshot is the restart fence used by fs-session:
    /// a fresh durable governor remains recovery-only until it has rebuilt the
    /// complete claim count observed at construction. Pending claims are
    /// included, so an indeterminate execution cannot be skipped to admit new
    /// work after a crash.
    ///
    /// # Errors
    /// Unavailable ledger identity, malformed aggregate storage, count
    /// overflow, or engine failure is returned without mutation.
    pub fn session_mutation_claim_count(
        &self,
        governor_hash: ContentHash,
    ) -> Result<u64, LedgerError> {
        let _ = self.checked_instance_id()?;
        self.note_read_query();
        let rows = self
            .conn
            .query_with_params(
                "SELECT \
                    (SELECT COUNT(*) FROM session_claims WHERE governor_hash = ?1), \
                    (SELECT COUNT(*) FROM session_claim_discovery WHERE governor_hash = ?1)",
                &[blob_param(governor_hash.as_bytes())],
            )
            .map_err(|error| sql_err("session governor claim count", &error))?;
        let row = rows.first().ok_or_else(|| {
            stored_corrupt(
                governor_hash,
                "governor claim-count query returned no aggregate row",
            )
        })?;
        let count = row_i64_registry(row, 0, governor_hash, "governor claim count")?;
        let discovery_count =
            row_i64_registry(row, 1, governor_hash, "governor claim discovery count")?;
        if count != discovery_count {
            return Err(stored_corrupt(
                governor_hash,
                format!(
                    "governor claim/discovery counts differ: claim={count}, discovery={discovery_count}"
                ),
            ));
        }
        u64::try_from(count).map_err(|_| {
            stored_corrupt(
                governor_hash,
                format!("governor claim count {count} is negative or overflows u64"),
            )
        })
    }

    /// Return one exact verified Pending claim in a recovery generation.
    ///
    /// The v7 composite index and bounded keyset pages make this a constrained
    /// recovery probe. A zero- or one-event singleton terminal is reauthenticated
    /// in one joined read; every other preceding terminalized claim is read
    /// through [`Self::session_terminal`]. Malformed raw state therefore cannot
    /// hide an indeterminate claim. The first fully verified Pending claim is
    /// returned.
    ///
    /// # Errors
    /// Invalid kind/scope bounds, malformed selected rows, unavailable ledger
    /// identity, or engine failures are returned without mutation.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)] // Complete bounded recovery-envelope verification is one probe.
    pub fn pending_session_mutation(
        &self,
        governor_hash: ContentHash,
        session_open_hash: ContentHash,
        kind: &str,
        session: u64,
        ledger_scope: &str,
        generation: u64,
    ) -> Result<Option<StoredSessionMutationClaim>, LedgerError> {
        require_bounded_ascii("session_claim.kind", kind, MAX_SESSION_TERMINAL_KIND_BYTES)?;
        require_bounded_ascii(
            "session_claim.ledger_scope",
            ledger_scope,
            MAX_SESSION_TERMINAL_SCOPE_BYTES,
        )?;
        let current_ledger = self.checked_instance_id()?;
        self.note_read_query();
        let simple_probe = self.prepare_simple_generation_probe()?;
        let mut cursor: Option<ContentHash> = None;
        let mut inspected = 0usize;
        loop {
            let mut params = vec![
                blob_param(governor_hash.as_bytes()),
                blob_param(session_open_hash.as_bytes()),
                text_param(kind),
                blob_param(&session.to_be_bytes()),
                text_param(ledger_scope),
                blob_param(&generation.to_be_bytes()),
            ];
            let sql = if let Some(after) = cursor {
                params.push(blob_param(after.as_bytes()));
                format!(
                    "SELECT authority FROM ( \
                         SELECT authority FROM session_claims \
                         WHERE governor_hash = ?1 AND session_open_hash = ?2 \
                           AND kind = ?3 AND session = ?4 AND ledger_scope = ?5 \
                           AND generation = ?6 AND authority > ?7 \
                         UNION \
                         SELECT authority FROM session_claim_discovery \
                         WHERE governor_hash = ?1 AND session_open_hash = ?2 \
                           AND kind = ?3 AND session = ?4 AND ledger_scope = ?5 \
                           AND generation = ?6 AND authority > ?7 \
                     ) AS discovered ORDER BY authority LIMIT {}",
                    MAX_SESSION_FLUSH_TERMINALS + 1
                )
            } else {
                format!(
                    "SELECT authority FROM ( \
                         SELECT authority FROM session_claims \
                         WHERE governor_hash = ?1 AND session_open_hash = ?2 \
                           AND kind = ?3 AND session = ?4 AND ledger_scope = ?5 \
                           AND generation = ?6 \
                         UNION \
                         SELECT authority FROM session_claim_discovery \
                         WHERE governor_hash = ?1 AND session_open_hash = ?2 \
                           AND kind = ?3 AND session = ?4 AND ledger_scope = ?5 \
                           AND generation = ?6 \
                     ) AS discovered ORDER BY authority LIMIT {}",
                    MAX_SESSION_FLUSH_TERMINALS + 1
                )
            };
            let rows = self
                .conn
                .query_with_params(&sql, &params)
                .map_err(|error| sql_err("session pending-claim probe", &error))?;
            let has_more = rows.len() > MAX_SESSION_FLUSH_TERMINALS;
            let mut last = None;
            for row in rows.iter().take(MAX_SESSION_FLUSH_TERMINALS) {
                inspected = inspected.checked_add(1).ok_or_else(|| {
                    invalid(
                        "session_claim.recovery_probe",
                        "generation recovery-probe count overflowed usize",
                    )
                })?;
                if inspected > MAX_SESSION_RECOVERY_PROBE_CLAIMS {
                    return Err(invalid(
                        "session_claim.recovery_probe",
                        format!(
                            "generation contains more than the {MAX_SESSION_RECOVERY_PROBE_CLAIMS}-claim recovery-probe limit"
                        ),
                    ));
                }
                let authority =
                    ContentHash(row_blob(row, 0, governor_hash, "pending claim authority")?);
                last = Some(authority);
                let (stored, terminalized) = if let Some(probe) =
                    self.simple_generation_probe_with_statement(&simple_probe, authority)?
                {
                    (probe.claim, Some(probe.terminalized))
                } else {
                    (
                        self.session_mutation_claim(&authority)?.ok_or_else(|| {
                            stored_corrupt(
                                authority,
                                "claim disappeared between indexed probe and verified read",
                            )
                        })?,
                        None,
                    )
                };
                if stored.ledger_instance_id != current_ledger
                    || stored.governor_hash != governor_hash
                    || stored.session_open_hash != session_open_hash
                    || stored.kind != kind
                    || stored.session != session
                    || stored.ledger_scope != ledger_scope
                    || stored.generation != generation
                {
                    return Err(stored_corrupt(
                        authority,
                        "claim differs from its indexed recovery envelope",
                    ));
                }
                let terminalized = match terminalized {
                    Some(terminalized) => terminalized,
                    None => self.session_terminal(&authority)?.is_some(),
                };
                if !terminalized {
                    return Ok(Some(stored));
                }
            }
            if !has_more {
                return Ok(None);
            }
            cursor = last;
        }
    }

    fn terminal_presence(&self, authority: &ContentHash) -> Result<bool, LedgerError> {
        let rows = self
            .conn
            .query_with_params(
                "SELECT 1 FROM session_terminals WHERE authority = ?1 LIMIT 1",
                &[blob_param(authority.as_bytes())],
            )
            .map_err(|error| sql_err("session terminal presence", &error))?;
        Ok(!rows.is_empty())
    }

    fn terminal_batch_memberships(
        &self,
        authority: &ContentHash,
    ) -> Result<Vec<(ContentHash, usize)>, LedgerError> {
        let sql = format!(
            "SELECT batch_id, seq FROM session_flush_batch_members \
             WHERE authority = ?1 ORDER BY batch_id LIMIT {}",
            MAX_SESSION_TERMINAL_BATCH_MEMBERSHIPS + 1
        );
        let rows = self
            .conn
            .query_with_params(&sql, &[blob_param(authority.as_bytes())])
            .map_err(|error| sql_err("session terminal batch memberships", &error))?;
        if rows.len() > MAX_SESSION_TERMINAL_BATCH_MEMBERSHIPS {
            return Err(stored_corrupt(
                *authority,
                format!(
                    "terminal authority exceeds the {MAX_SESSION_TERMINAL_BATCH_MEMBERSHIPS}-batch witness limit"
                ),
            ));
        }
        rows.iter()
            .map(|row| {
                Ok::<_, LedgerError>((
                    ContentHash(row_blob(row, 0, *authority, "batch_member.batch_id")?),
                    row_usize(row, 1, *authority, "batch_member.seq")?,
                ))
            })
            .collect()
    }

    fn terminal_event_link_count(&self, authority: ContentHash) -> Result<usize, LedgerError> {
        let rows = self
            .conn
            .query_with_params(
                "SELECT COUNT(*) FROM session_terminal_events WHERE authority = ?1",
                &[blob_param(authority.as_bytes())],
            )
            .map_err(|error| sql_err("session terminal event-link count", &error))?;
        let row = rows
            .first()
            .ok_or_else(|| stored_corrupt(authority, "event-link count query returned no row"))?;
        row_usize(row, 0, authority, "event_link_count")
    }

    #[allow(clippy::too_many_lines)] // One bounded metadata/read/hash verification protocol.
    fn verify_session_terminal_events(
        &self,
        claim: &StoredSessionMutationClaim,
        receipt_len: usize,
        expected_count: usize,
        expected_hash: ContentHash,
        expected_encoded_bytes: usize,
    ) -> Result<(), LedgerError> {
        let authority = claim.authority;
        let link_count = self.terminal_event_link_count(authority)?;
        if link_count != expected_count || link_count > MAX_SESSION_FLUSH_EVENTS {
            return Err(stored_corrupt(
                authority,
                format!(
                    "terminal records {expected_count} events but has {link_count} authority links"
                ),
            ));
        }

        let metadata_sql = format!(
            "SELECT link.seq, link.event_id, \
                    CASE WHEN owned.id IS NULL THEN 0 ELSE 1 END, \
                    CASE WHEN typeof(owned.session) = 'blob' \
                         THEN length(owned.session) ELSE -1 END, \
                    CASE WHEN typeof(owned.t) = 'integer' THEN 1 ELSE 0 END, \
                    CASE WHEN typeof(owned.kind) = 'text' AND \
                                   length(CAST(owned.kind AS BLOB)) BETWEEN 1 AND \
                                       {MAX_SESSION_TERMINAL_EVENT_KIND_BYTES} AND \
                                   length(CAST(owned.kind AS BLOB)) = length(owned.kind) AND \
                                   owned.kind NOT GLOB '*[^!-~]*' \
                         THEN length(CAST(owned.kind AS BLOB)) ELSE -1 END, \
                    CASE WHEN owned.payload IS NULL THEN 0 \
                         WHEN typeof(owned.payload) = 'text' AND \
                              length(CAST(owned.payload AS BLOB)) BETWEEN 1 AND \
                                  {MAX_SESSION_TERMINAL_EVENT_PAYLOAD_BYTES} \
                         THEN length(CAST(owned.payload AS BLOB)) + 1 ELSE -1 END \
             FROM session_terminal_events AS link \
             LEFT JOIN events AS owned ON owned.id = link.event_id \
             WHERE link.authority = ?1 ORDER BY link.seq LIMIT {}",
            MAX_SESSION_FLUSH_EVENTS + 1
        );
        let metadata = self
            .conn
            .query_with_params(&metadata_sql, &[blob_param(authority.as_bytes())])
            .map_err(|error| sql_err("session terminal event metadata", &error))?;
        if metadata.len() != expected_count {
            return Err(stored_corrupt(
                authority,
                "event-link metadata count changed during bounded verification",
            ));
        }

        let mut recomputed_encoded_bytes = SESSION_CLAIM_ROW_FRAMING_BYTES
            .checked_add(claim.kind.len())
            .and_then(|bytes| bytes.checked_add(claim.ledger_scope.len()))
            .and_then(|bytes| bytes.checked_add(claim.payload.len()))
            .and_then(|bytes| bytes.checked_add(SESSION_TERMINAL_ROW_FRAMING_BYTES))
            .and_then(|bytes| bytes.checked_add(receipt_len))
            .ok_or_else(|| {
                stored_corrupt(
                    authority,
                    "stored terminal encoded-byte count overflowed usize",
                )
            })?;
        for (expected_seq, row) in metadata.iter().enumerate() {
            let seq = row_usize(row, 0, authority, "event_link.seq")?;
            let event_id = row_i64_registry(row, 1, authority, "event_link.event_id")?;
            let present = row_i64_registry(row, 2, authority, "event presence")?;
            let session_len = row_i64_registry(row, 3, authority, "event session length")?;
            let t_is_integer = row_i64_registry(row, 4, authority, "event t type")?;
            let kind_len = row_i64_registry(row, 5, authority, "event kind length")?;
            let payload_state = row_i64_registry(row, 6, authority, "event payload state")?;
            if seq != expected_seq
                || event_id <= 0
                || present != 1
                || session_len != 8
                || t_is_integer != 1
                || kind_len <= 0
                || payload_state < 0
            {
                return Err(stored_corrupt(
                    authority,
                    format!("event ownership metadata is malformed at sequence {expected_seq}"),
                ));
            }
            let payload_len = if payload_state == 0 {
                0
            } else {
                usize::try_from(payload_state - 1).map_err(|_| {
                    stored_corrupt(authority, "event payload length does not fit usize")
                })?
            };
            let kind_len = usize::try_from(kind_len)
                .map_err(|_| stored_corrupt(authority, "event kind length does not fit usize"))?;
            let event_bytes = SESSION_EVENT_ROW_FRAMING_BYTES
                .checked_add(8)
                .and_then(|bytes| bytes.checked_add(kind_len))
                .and_then(|bytes| bytes.checked_add(payload_len))
                .ok_or_else(|| {
                    stored_corrupt(authority, "event encoded-byte count overflowed usize")
                })?;
            recomputed_encoded_bytes = recomputed_encoded_bytes
                .checked_add(event_bytes)
                .ok_or_else(|| {
                    stored_corrupt(authority, "terminal encoded-byte count overflowed usize")
                })?;
            if recomputed_encoded_bytes > MAX_SESSION_FLUSH_ENCODED_BYTES {
                return Err(stored_corrupt(
                    authority,
                    "terminal exceeds the aggregate encoded-byte bound",
                ));
            }
        }
        if recomputed_encoded_bytes != expected_encoded_bytes {
            return Err(stored_corrupt(
                authority,
                format!(
                    "terminal encoded-byte count is {expected_encoded_bytes}, recomputed \
                     {recomputed_encoded_bytes}"
                ),
            ));
        }

        let guarded_sql = format!(
            "SELECT link.seq, owned.session, owned.t, owned.kind, owned.payload \
             FROM session_terminal_events AS link \
             JOIN events AS owned ON owned.id = link.event_id \
             WHERE link.authority = ?1 AND \
               typeof(link.seq) = 'integer' AND \
               typeof(link.event_id) = 'integer' AND \
               typeof(owned.session) = 'blob' AND length(owned.session) = 8 AND \
               typeof(owned.t) = 'integer' AND \
               typeof(owned.kind) = 'text' AND \
               length(CAST(owned.kind AS BLOB)) BETWEEN 1 AND \
                   {MAX_SESSION_TERMINAL_EVENT_KIND_BYTES} AND \
               length(CAST(owned.kind AS BLOB)) = length(owned.kind) AND \
               owned.kind NOT GLOB '*[^!-~]*' AND \
               CASE WHEN owned.payload IS NULL THEN 1 \
                    WHEN typeof(owned.payload) = 'text' AND \
                         length(CAST(owned.payload AS BLOB)) BETWEEN 1 AND \
                             {MAX_SESSION_TERMINAL_EVENT_PAYLOAD_BYTES} \
                    THEN json_valid(owned.payload) ELSE 0 END = 1 \
             ORDER BY link.seq LIMIT {}",
            MAX_SESSION_FLUSH_EVENTS + 1
        );
        let rows = self
            .conn
            .query_with_params(&guarded_sql, &[blob_param(authority.as_bytes())])
            .map_err(|error| sql_err("session terminal owned-event fetch", &error))?;
        if rows.len() != expected_count {
            return Err(stored_corrupt(
                authority,
                "owned global event is missing or violates a type, JSON, or byte bound",
            ));
        }

        let expected_session = claim.session.to_be_bytes();
        let mut hasher = events_hasher(expected_count);
        for (expected_seq, row) in rows.iter().enumerate() {
            if row_usize(row, 0, authority, "event_link.seq")? != expected_seq {
                return Err(stored_corrupt(
                    authority,
                    format!("event link sequence is not dense at {expected_seq}"),
                ));
            }
            let session: [u8; 8] = row_blob(row, 1, authority, "event.session")?;
            if session != expected_session {
                return Err(stored_corrupt(
                    authority,
                    format!("owned event {expected_seq} belongs to another session"),
                ));
            }
            let t = row_i64_registry(row, 2, authority, "event.t")?;
            let Some(SqliteValue::Text(kind)) = row.get(3) else {
                return Err(stored_corrupt(authority, "event.kind is not TEXT"));
            };
            let payload = match row.get(4) {
                Some(SqliteValue::Null) => None,
                Some(SqliteValue::Text(payload)) => Some(payload.as_str()),
                _ => {
                    return Err(stored_corrupt(
                        authority,
                        "event.payload is not NULL or TEXT",
                    ));
                }
            };
            update_event_preimage(&mut hasher, &session, t, kind.as_str(), payload);
        }
        if hasher.finalize() != expected_hash {
            return Err(stored_corrupt(
                authority,
                "owned global event bytes do not match events_hash",
            ));
        }
        Ok(())
    }

    /// Fetch one terminal receipt and verify its claim, receipt hash, dense
    /// ownership links, joined global events, aggregate bytes, and event hash.
    ///
    /// Absence is Ok(None) for an absent or still-Pending claim.
    ///
    /// # Errors
    /// Any orphan, partial, foreign-ledger, malformed, or hash-mismatched state
    /// fails closed.
    pub fn session_terminal(
        &self,
        authority: &ContentHash,
    ) -> Result<Option<StoredSessionTerminal>, LedgerError> {
        self.note_read_query();
        let Some(claim) = self.session_mutation_claim(authority)? else {
            let terminal_present = self.terminal_presence(authority)?;
            let batch_memberships = self.terminal_batch_memberships(authority)?;
            let link_count = self.terminal_event_link_count(*authority)?;
            return if terminal_present || !batch_memberships.is_empty() || link_count != 0 {
                Err(stored_corrupt(
                    *authority,
                    "terminal ownership exists without its immutable mutation claim",
                ))
            } else {
                Ok(None)
            };
        };
        let terminal_present = self.terminal_presence(authority)?;
        let batch_memberships = self.terminal_batch_memberships(authority)?;
        if !terminal_present {
            let link_count = self.terminal_event_link_count(*authority)?;
            return if !batch_memberships.is_empty() || link_count != 0 {
                Err(stored_corrupt(
                    *authority,
                    "terminal row is missing but durable batch or event ownership remains",
                ))
            } else {
                Ok(None)
            };
        }
        if batch_memberships.is_empty() {
            return Err(stored_corrupt(
                *authority,
                "terminal row has no immutable flush-batch membership witness",
            ));
        }

        let guarded_sql = format!(
            "SELECT receipt, receipt_hash, event_count, events_hash, encoded_bytes, created_at \
             FROM session_terminals WHERE authority = ?1 AND \
               typeof(receipt) = 'blob' AND \
               length(receipt) BETWEEN 1 AND {MAX_SESSION_TERMINAL_RECEIPT_BYTES} AND \
               typeof(receipt_hash) = 'blob' AND length(receipt_hash) = 32 AND \
               typeof(event_count) = 'integer' AND \
               event_count BETWEEN 0 AND {MAX_SESSION_FLUSH_EVENTS} AND \
               typeof(events_hash) = 'blob' AND length(events_hash) = 32 AND \
               typeof(encoded_bytes) = 'integer' AND \
               encoded_bytes BETWEEN 1 AND {MAX_SESSION_FLUSH_ENCODED_BYTES} AND \
               typeof(created_at) = 'integer' LIMIT 1"
        );
        let rows = self
            .conn
            .query_with_params(&guarded_sql, &[blob_param(authority.as_bytes())])
            .map_err(|error| sql_err("session terminal guarded fetch", &error))?;
        let row = rows.first().ok_or_else(|| {
            stored_corrupt(
                *authority,
                "terminal violates a type, hash, count, or byte bound",
            )
        })?;
        let receipt = row_blob_vec(row, 0, *authority, "receipt")?;
        let receipt_hash = ContentHash(row_blob(row, 1, *authority, "receipt_hash")?);
        if hash_bytes(&receipt) != receipt_hash {
            return Err(stored_corrupt(
                *authority,
                "receipt bytes do not match receipt_hash",
            ));
        }
        let event_count = row_usize(row, 2, *authority, "event_count")?;
        let events_hash = ContentHash(row_blob(row, 3, *authority, "events_hash")?);
        let encoded_bytes = row_usize(row, 4, *authority, "encoded_bytes")?;
        self.verify_session_terminal_events(
            &claim,
            receipt.len(),
            event_count,
            events_hash,
            encoded_bytes,
        )?;
        for (batch_id, batch_seq) in batch_memberships {
            let batch = self.session_flush_batch(&batch_id)?.ok_or_else(|| {
                stored_corrupt(
                    *authority,
                    format!("terminal membership references missing flush batch {batch_id}"),
                )
            })?;
            if batch_seq >= batch.terminal_count {
                return Err(stored_corrupt(
                    *authority,
                    format!(
                        "terminal membership sequence {batch_seq} exceeds batch {} terminal count {}",
                        batch.batch_id, batch.terminal_count
                    ),
                ));
            }
        }
        Ok(Some(StoredSessionTerminal {
            claim,
            receipt,
            receipt_hash,
            event_count,
            events_hash,
            encoded_bytes,
            created_at: row_i64_registry(row, 5, *authority, "created_at")?,
        }))
    }

    fn verify_session_flush_batch_members(
        &self,
        batch: &StoredSessionFlushBatch,
    ) -> Result<(), LedgerError> {
        let sql = format!(
            "SELECT member.seq, member.authority, \
                    CASE WHEN terminal.authority IS NOT NULL AND claim.authority IS NOT NULL \
                         THEN 1 ELSE 0 END, \
                    claim.claim_hash, terminal.receipt_hash, terminal.event_count, \
                    terminal.events_hash, terminal.encoded_bytes \
             FROM session_flush_batch_members AS member \
             LEFT JOIN session_terminals AS terminal \
                    ON terminal.authority = member.authority \
             LEFT JOIN session_claims AS claim ON claim.authority = member.authority \
             WHERE member.batch_id = ?1 ORDER BY member.seq LIMIT {}",
            MAX_SESSION_FLUSH_TERMINALS + 1
        );
        let rows = self
            .conn
            .query_with_params(&sql, &[blob_param(batch.batch_id.as_bytes())])
            .map_err(|error| sql_err("session flush batch-member verification", &error))?;
        if rows.len() != batch.terminal_count {
            return Err(stored_corrupt(
                batch.batch_id,
                format!(
                    "batch records {} terminals but has {} ordered membership rows",
                    batch.terminal_count,
                    rows.len()
                ),
            ));
        }
        let mut authorities = BTreeSet::new();
        let mut event_count = 0usize;
        let mut encoded_bytes = 0usize;
        let mut hasher = Blake3::new();
        hasher.update(SESSION_BATCH_HASH_DOMAIN);
        hasher.update(&batch.ledger_instance_id.as_bytes());
        hasher.update(&SESSION_REGISTRY_ROW_SCHEMA_VERSION.to_le_bytes());
        update_len(&mut hasher, batch.terminal_count);
        for (expected_seq, row) in rows.iter().enumerate() {
            if row_usize(row, 0, batch.batch_id, "batch_member.seq")? != expected_seq {
                return Err(stored_corrupt(
                    batch.batch_id,
                    format!("batch membership sequence is not dense at {expected_seq}"),
                ));
            }
            let authority =
                ContentHash(row_blob(row, 1, batch.batch_id, "batch_member.authority")?);
            if !authorities.insert(authority) {
                return Err(stored_corrupt(
                    batch.batch_id,
                    format!("authority {authority} appears twice in one batch"),
                ));
            }
            if row_i64_registry(row, 2, batch.batch_id, "batch member target presence")? != 1 {
                return Err(stored_corrupt(
                    batch.batch_id,
                    format!("batch member {expected_seq} lacks its claim or terminal row"),
                ));
            }
            let claim_hash =
                ContentHash(row_blob(row, 3, batch.batch_id, "batch_member.claim_hash")?);
            let receipt_hash = ContentHash(row_blob(
                row,
                4,
                batch.batch_id,
                "batch_member.receipt_hash",
            )?);
            let terminal_events = row_usize(row, 5, batch.batch_id, "batch_member.event_count")?;
            if terminal_events > MAX_SESSION_FLUSH_EVENTS {
                return Err(stored_corrupt(
                    batch.batch_id,
                    format!(
                        "batch member {expected_seq} records {terminal_events} events above the {MAX_SESSION_FLUSH_EVENTS} limit"
                    ),
                ));
            }
            let events_hash = ContentHash(row_blob(
                row,
                6,
                batch.batch_id,
                "batch_member.events_hash",
            )?);
            let terminal_bytes = row_usize(row, 7, batch.batch_id, "batch_member.encoded_bytes")?;
            if terminal_bytes == 0 || terminal_bytes > MAX_SESSION_FLUSH_ENCODED_BYTES {
                return Err(stored_corrupt(
                    batch.batch_id,
                    format!(
                        "batch member {expected_seq} encoded byte count {terminal_bytes} is outside 1..={MAX_SESSION_FLUSH_ENCODED_BYTES}"
                    ),
                ));
            }
            event_count = event_count.checked_add(terminal_events).ok_or_else(|| {
                stored_corrupt(
                    batch.batch_id,
                    "batch member event-count sum overflowed usize",
                )
            })?;
            encoded_bytes = encoded_bytes.checked_add(terminal_bytes).ok_or_else(|| {
                stored_corrupt(
                    batch.batch_id,
                    "batch member byte-count sum overflowed usize",
                )
            })?;
            hasher.update(authority.as_bytes());
            hasher.update(claim_hash.as_bytes());
            hasher.update(receipt_hash.as_bytes());
            update_len(&mut hasher, terminal_events);
            hasher.update(events_hash.as_bytes());
            update_len(&mut hasher, terminal_bytes);
        }
        if event_count != batch.event_count || encoded_bytes != batch.encoded_bytes {
            return Err(stored_corrupt(
                batch.batch_id,
                format!(
                    "batch marker totals ({}, {}) disagree with member totals ({event_count}, {encoded_bytes})",
                    batch.event_count, batch.encoded_bytes
                ),
            ));
        }
        if hasher.finalize() != batch.batch_id {
            return Err(stored_corrupt(
                batch.batch_id,
                "batch id does not authenticate its complete ordered membership preimage",
            ));
        }
        Ok(())
    }

    /// Fetch one immutable canonical flush-batch marker.
    ///
    /// # Errors
    /// Foreign-ledger, future-envelope, malformed, or out-of-bound rows fail
    /// closed.
    pub fn session_flush_batch(
        &self,
        batch_id: &ContentHash,
    ) -> Result<Option<StoredSessionFlushBatch>, LedgerError> {
        let current_ledger = self.checked_instance_id()?;
        self.note_read_query();
        let present = self
            .conn
            .query_with_params(
                "SELECT 1 FROM session_flush_batches WHERE batch_id = ?1 LIMIT 1",
                &[blob_param(batch_id.as_bytes())],
            )
            .map_err(|error| sql_err("session flush batch presence", &error))?;
        if present.is_empty() {
            return Ok(None);
        }
        let guarded_sql = format!(
            "SELECT ledger_instance_id, registry_schema_version, terminal_count, event_count, \
                    encoded_bytes, created_at \
             FROM session_flush_batches WHERE batch_id = ?1 AND \
               typeof(ledger_instance_id) = 'blob' AND length(ledger_instance_id) = 16 AND \
               typeof(registry_schema_version) = 'integer' AND \
               registry_schema_version = {SESSION_REGISTRY_ROW_SCHEMA_VERSION} AND \
               typeof(terminal_count) = 'integer' AND \
               terminal_count BETWEEN 1 AND {MAX_SESSION_FLUSH_TERMINALS} AND \
               typeof(event_count) = 'integer' AND \
               event_count BETWEEN 0 AND {MAX_SESSION_FLUSH_EVENTS} AND \
               typeof(encoded_bytes) = 'integer' AND \
               encoded_bytes BETWEEN 1 AND {MAX_SESSION_FLUSH_ENCODED_BYTES} AND \
               typeof(created_at) = 'integer' LIMIT 1"
        );
        let rows = self
            .conn
            .query_with_params(&guarded_sql, &[blob_param(batch_id.as_bytes())])
            .map_err(|error| sql_err("session flush batch guarded fetch", &error))?;
        let row = rows.first().ok_or_else(|| {
            stored_corrupt(
                *batch_id,
                "flush-batch row violates a type, schema-version, or numeric bound",
            )
        })?;
        let stored_ledger: [u8; 16] = row_blob(row, 0, *batch_id, "ledger_instance_id")?;
        if stored_ledger != current_ledger.as_bytes() {
            return Err(stored_corrupt(
                *batch_id,
                "flush-batch marker belongs to a different physical ledger",
            ));
        }
        let stored = StoredSessionFlushBatch {
            batch_id: *batch_id,
            ledger_instance_id: current_ledger,
            schema_version: row_i64_registry(row, 1, *batch_id, "registry_schema_version")?,
            terminal_count: row_usize(row, 2, *batch_id, "terminal_count")?,
            event_count: row_usize(row, 3, *batch_id, "event_count")?,
            encoded_bytes: row_usize(row, 4, *batch_id, "encoded_bytes")?,
            created_at: row_i64_registry(row, 5, *batch_id, "created_at")?,
        };
        self.verify_session_flush_batch_members(&stored)?;
        Ok(Some(stored))
    }

    fn insert_session_claim(
        &self,
        claim: SessionMutationClaim<'_>,
        prepared: PreparedClaim,
    ) -> Result<(), LedgerError> {
        self.conn
            .prepare(
                "INSERT INTO session_claims( \
                    authority, ledger_instance_id, governor_hash, session_open_hash, \
                    registry_schema_version, kind, session, ledger_scope, generation, \
                    causal_ordinal, payload, payload_hash, claim_hash, created_at \
                 ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13, ?14)",
            )
            .map_err(|error| sql_err("session claim insert prepare", &error))?
            .execute_with_params(&[
                blob_param(claim.authority.as_bytes()),
                blob_param(&claim.ledger_instance_id.as_bytes()),
                blob_param(claim.governor_hash.as_bytes()),
                blob_param(claim.session_open_hash.as_bytes()),
                SqliteValue::Integer(SESSION_REGISTRY_ROW_SCHEMA_VERSION),
                text_param(claim.kind),
                blob_param(&claim.session.to_be_bytes()),
                text_param(claim.ledger_scope),
                blob_param(&claim.generation.to_be_bytes()),
                claim.causal_ordinal.map_or(SqliteValue::Null, |ordinal| {
                    blob_param(&ordinal.to_be_bytes())
                }),
                blob_param(claim.payload),
                blob_param(prepared.payload_hash.as_bytes()),
                blob_param(prepared.claim_hash.as_bytes()),
                SqliteValue::Integer(now_wall_ns()),
            ])
            .map_err(|error| sql_err("session claim insert", &error))?;
        self.conn
            .prepare(
                "INSERT INTO session_claim_discovery( \
                    authority, ledger_instance_id, governor_hash, session_open_hash, \
                    registry_schema_version, kind, session, ledger_scope, generation, \
                    causal_ordinal, payload_hash, claim_hash \
                 ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12)",
            )
            .map_err(|error| sql_err("session claim discovery insert prepare", &error))?
            .execute_with_params(&[
                blob_param(claim.authority.as_bytes()),
                blob_param(&claim.ledger_instance_id.as_bytes()),
                blob_param(claim.governor_hash.as_bytes()),
                blob_param(claim.session_open_hash.as_bytes()),
                SqliteValue::Integer(SESSION_REGISTRY_ROW_SCHEMA_VERSION),
                text_param(claim.kind),
                blob_param(&claim.session.to_be_bytes()),
                text_param(claim.ledger_scope),
                blob_param(&claim.generation.to_be_bytes()),
                claim.causal_ordinal.map_or(SqliteValue::Null, |ordinal| {
                    blob_param(&ordinal.to_be_bytes())
                }),
                blob_param(prepared.payload_hash.as_bytes()),
                blob_param(prepared.claim_hash.as_bytes()),
            ])
            .map_err(|error| sql_err("session claim discovery insert", &error))?;
        Ok(())
    }

    fn insert_session_terminal(
        &self,
        row: SessionTerminalRow<'_>,
        prepared: PreparedTerminal,
    ) -> Result<(), LedgerError> {
        self.conn
            .prepare(
                "INSERT INTO session_terminals( \
                    authority, receipt, receipt_hash, event_count, events_hash, \
                    encoded_bytes, created_at \
                 ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)",
            )
            .map_err(|error| sql_err("session terminal insert prepare", &error))?
            .execute_with_params(&[
                blob_param(row.claim.authority.as_bytes()),
                blob_param(row.receipt),
                blob_param(prepared.receipt_hash.as_bytes()),
                SqliteValue::Integer(
                    i64::try_from(prepared.event_count)
                        .expect("bounded terminal event count fits i64"),
                ),
                blob_param(prepared.events_hash.as_bytes()),
                SqliteValue::Integer(
                    i64::try_from(prepared.encoded_bytes)
                        .expect("bounded terminal encoded bytes fit i64"),
                ),
                SqliteValue::Integer(now_wall_ns()),
            ])
            .map_err(|error| sql_err("session terminal insert", &error))?;
        Ok(())
    }

    fn insert_session_terminal_event_link(
        &self,
        authority: ContentHash,
        seq: usize,
        event: i64,
    ) -> Result<(), LedgerError> {
        self.conn
            .prepare(
                "INSERT INTO session_terminal_events(authority, seq, event_id) \
                 VALUES (?1, ?2, ?3)",
            )
            .map_err(|error| sql_err("session terminal event-link insert prepare", &error))?
            .execute_with_params(&[
                blob_param(authority.as_bytes()),
                SqliteValue::Integer(
                    i64::try_from(seq).expect("bounded terminal event sequence fits i64"),
                ),
                SqliteValue::Integer(event),
            ])
            .map_err(|error| sql_err("session terminal event-link insert", &error))?;
        Ok(())
    }

    fn insert_session_flush_batch(
        &self,
        prepared: &PreparedBatch,
        terminal_count: usize,
    ) -> Result<(), LedgerError> {
        self.conn
            .prepare(
                "INSERT INTO session_flush_batches( \
                    batch_id, ledger_instance_id, registry_schema_version, terminal_count, \
                    event_count, encoded_bytes, created_at \
                 ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)",
            )
            .map_err(|error| sql_err("session flush batch insert prepare", &error))?
            .execute_with_params(&[
                blob_param(prepared.batch_id.as_bytes()),
                blob_param(&prepared.ledger_instance_id.as_bytes()),
                SqliteValue::Integer(SESSION_REGISTRY_ROW_SCHEMA_VERSION),
                SqliteValue::Integer(
                    i64::try_from(terminal_count).expect("bounded terminal count fits i64"),
                ),
                SqliteValue::Integer(
                    i64::try_from(prepared.event_count).expect("bounded event count fits i64"),
                ),
                SqliteValue::Integer(
                    i64::try_from(prepared.encoded_bytes)
                        .expect("bounded encoded byte count fits i64"),
                ),
                SqliteValue::Integer(now_wall_ns()),
            ])
            .map_err(|error| sql_err("session flush batch insert", &error))?;
        Ok(())
    }

    fn insert_session_flush_batch_member(
        &self,
        batch_id: ContentHash,
        seq: usize,
        authority: ContentHash,
    ) -> Result<(), LedgerError> {
        self.conn
            .prepare(
                "INSERT INTO session_flush_batch_members(batch_id, seq, authority) \
                 VALUES (?1, ?2, ?3)",
            )
            .map_err(|error| sql_err("session flush batch-member insert prepare", &error))?
            .execute_with_params(&[
                blob_param(batch_id.as_bytes()),
                SqliteValue::Integer(
                    i64::try_from(seq).expect("bounded batch-member sequence fits i64"),
                ),
                blob_param(authority.as_bytes()),
            ])
            .map_err(|error| sql_err("session flush batch-member insert", &error))?;
        Ok(())
    }

    fn terminalized_pause_successor(
        &self,
        claim: SessionMutationClaim<'_>,
    ) -> Result<Option<ContentHash>, LedgerError> {
        let Some(resume_generation) = claim.generation.checked_add(1) else {
            return Ok(None);
        };
        let rows = self
            .conn
            .query_with_params(
                "SELECT pause.authority FROM ( \
                     SELECT authority FROM session_claims \
                     WHERE governor_hash = ?1 AND session_open_hash = ?2 \
                       AND kind = ?3 AND session = ?4 \
                       AND ledger_scope = ?5 AND generation = ?6 \
                     UNION \
                     SELECT authority FROM session_claim_discovery \
                     WHERE governor_hash = ?1 AND session_open_hash = ?2 \
                       AND kind = ?3 AND session = ?4 \
                       AND ledger_scope = ?5 AND generation = ?6 \
                 ) AS pause WHERE ( \
                       EXISTS(SELECT 1 FROM session_terminals AS terminal \
                              WHERE terminal.authority = pause.authority) \
                       OR EXISTS(SELECT 1 FROM session_flush_batch_members AS member \
                                 WHERE member.authority = pause.authority) \
                       OR EXISTS(SELECT 1 FROM session_terminal_events AS event_link \
                                 WHERE event_link.authority = pause.authority) \
                   ) \
                 ORDER BY pause.authority LIMIT 1",
                &[
                    blob_param(claim.governor_hash.as_bytes()),
                    blob_param(claim.session_open_hash.as_bytes()),
                    text_param(PAUSE_ACKNOWLEDGEMENT_KIND),
                    blob_param(&claim.session.to_be_bytes()),
                    text_param(claim.ledger_scope),
                    blob_param(&resume_generation.to_be_bytes()),
                ],
            )
            .map_err(|error| sql_err("session successor pause fence", &error))?;
        let Some(row) = rows.first() else {
            return Ok(None);
        };
        let authority = ContentHash(row_blob(
            row,
            0,
            claim.authority,
            "successor pause authority",
        )?);
        let stored = self.session_mutation_claim(&authority)?.ok_or_else(|| {
            stored_corrupt(
                authority,
                "successor pause claim disappeared during generation-fence verification",
            )
        })?;
        if stored.governor_hash != claim.governor_hash
            || stored.session_open_hash != claim.session_open_hash
            || stored.kind != PAUSE_ACKNOWLEDGEMENT_KIND
            || stored.session != claim.session
            || stored.ledger_scope != claim.ledger_scope
            || stored.generation != resume_generation
        {
            return Err(stored_corrupt(
                authority,
                "successor pause claim differs from its generation-fence index envelope",
            ));
        }
        Ok(self.session_terminal(&authority)?.map(|_| authority))
    }

    fn causal_ordinal_owner(
        &self,
        claim: SessionMutationClaim<'_>,
    ) -> Result<Option<ContentHash>, LedgerError> {
        let Some(ordinal) = claim.causal_ordinal else {
            return Ok(None);
        };
        let rows = self
            .conn
            .query_with_params(
                "SELECT authority FROM ( \
                     SELECT authority FROM session_claims \
                     WHERE governor_hash = ?1 AND kind = ?2 AND causal_ordinal = ?3 \
                     UNION \
                     SELECT authority FROM session_claim_discovery \
                     WHERE governor_hash = ?1 AND kind = ?2 AND causal_ordinal = ?3 \
                 ) AS discovered ORDER BY authority LIMIT 2",
                &[
                    blob_param(claim.governor_hash.as_bytes()),
                    text_param(claim.kind),
                    blob_param(&ordinal.to_be_bytes()),
                ],
            )
            .map_err(|error| sql_err("session causal-ordinal owner", &error))?;
        if rows.len() > 1 {
            return Err(stored_corrupt(
                claim.authority,
                format!("governor/kind causal ordinal {ordinal} has multiple immutable owners"),
            ));
        }
        let Some(row) = rows.first() else {
            return Ok(None);
        };
        let authority = ContentHash(row_blob(row, 0, claim.authority, "causal ordinal owner")?);
        let stored = self.session_mutation_claim(&authority)?.ok_or_else(|| {
            stored_corrupt(
                authority,
                "causal-ordinal discovery witness references a missing claim",
            )
        })?;
        if stored.governor_hash != claim.governor_hash
            || stored.kind != claim.kind
            || stored.causal_ordinal != Some(ordinal)
        {
            return Err(stored_corrupt(
                authority,
                format!("causal-ordinal discovery differs from governor/kind ordinal {ordinal}"),
            ));
        }
        Ok(Some(authority))
    }

    #[allow(clippy::too_many_lines)] // Bounded keyset scan verifies every predecessor before one atomic fence decision.
    fn pending_submission_predecessors(
        &self,
        pause: SessionMutationClaim<'_>,
    ) -> Result<BTreeSet<ContentHash>, LedgerError> {
        let Some(draining_generation) = pause.generation.checked_sub(1) else {
            return Ok(BTreeSet::new());
        };
        let simple_probe = self.prepare_simple_generation_probe()?;
        let mut authorities = BTreeSet::new();
        let mut cursor: Option<ContentHash> = None;
        let mut inspected = 0usize;
        loop {
            let mut params = vec![
                blob_param(pause.governor_hash.as_bytes()),
                blob_param(pause.session_open_hash.as_bytes()),
                text_param(PRECLAIM_REQUIRED_SUBMISSION_KIND),
                blob_param(&pause.session.to_be_bytes()),
                text_param(pause.ledger_scope),
                blob_param(&draining_generation.to_be_bytes()),
            ];
            let sql = if let Some(after) = cursor {
                params.push(blob_param(after.as_bytes()));
                format!(
                    "SELECT authority FROM ( \
                         SELECT authority FROM session_claims \
                         WHERE governor_hash = ?1 AND session_open_hash = ?2 \
                           AND kind = ?3 AND session = ?4 AND ledger_scope = ?5 \
                           AND generation = ?6 AND authority > ?7 \
                         UNION \
                         SELECT authority FROM session_claim_discovery \
                         WHERE governor_hash = ?1 AND session_open_hash = ?2 \
                           AND kind = ?3 AND session = ?4 AND ledger_scope = ?5 \
                           AND generation = ?6 AND authority > ?7 \
                     ) AS discovered ORDER BY authority LIMIT {}",
                    MAX_SESSION_FLUSH_TERMINALS + 1
                )
            } else {
                format!(
                    "SELECT authority FROM ( \
                         SELECT authority FROM session_claims \
                         WHERE governor_hash = ?1 AND session_open_hash = ?2 \
                           AND kind = ?3 AND session = ?4 AND ledger_scope = ?5 \
                           AND generation = ?6 \
                         UNION \
                         SELECT authority FROM session_claim_discovery \
                         WHERE governor_hash = ?1 AND session_open_hash = ?2 \
                           AND kind = ?3 AND session = ?4 AND ledger_scope = ?5 \
                           AND generation = ?6 \
                     ) AS discovered ORDER BY authority LIMIT {}",
                    MAX_SESSION_FLUSH_TERMINALS + 1
                )
            };
            let rows = self
                .conn
                .query_with_params(&sql, &params)
                .map_err(|error| sql_err("session pending-submission pause fence", &error))?;
            let has_more = rows.len() > MAX_SESSION_FLUSH_TERMINALS;
            let page = rows.iter().take(MAX_SESSION_FLUSH_TERMINALS);
            let mut last = None;
            for row in page {
                inspected = inspected.checked_add(1).ok_or_else(|| {
                    invalid(
                        "session_terminal.pause_fence",
                        "submission predecessor count overflowed usize",
                    )
                })?;
                if inspected > MAX_SESSION_PAUSE_FENCE_SUBMISSIONS {
                    return Err(invalid(
                        "session_terminal.pause_fence",
                        format!(
                            "draining generation contains more than the {MAX_SESSION_PAUSE_FENCE_SUBMISSIONS}-submission fence limit"
                        ),
                    ));
                }
                let authority = ContentHash(row_blob(
                    row,
                    0,
                    pause.authority,
                    "submission predecessor authority",
                )?);
                last = Some(authority);
                let (stored, terminalized) = if let Some(probe) =
                    self.simple_generation_probe_with_statement(&simple_probe, authority)?
                {
                    (probe.claim, Some(probe.terminalized))
                } else {
                    (
                            self.session_mutation_claim(&authority)?.ok_or_else(|| {
                                stored_corrupt(
                                    authority,
                                    "submission predecessor disappeared during pause-fence verification",
                                )
                            })?,
                            None,
                        )
                };
                if stored.governor_hash != pause.governor_hash
                    || stored.session_open_hash != pause.session_open_hash
                    || stored.kind != PRECLAIM_REQUIRED_SUBMISSION_KIND
                    || stored.session != pause.session
                    || stored.ledger_scope != pause.ledger_scope
                    || stored.generation != draining_generation
                {
                    return Err(stored_corrupt(
                        authority,
                        "submission predecessor differs from its pause-fence index envelope",
                    ));
                }
                let terminalized = match terminalized {
                    Some(terminalized) => terminalized,
                    None => self.session_terminal(&authority)?.is_some(),
                };
                if !terminalized {
                    authorities.insert(authority);
                    if authorities.len() > MAX_SESSION_FLUSH_TERMINALS {
                        return Err(invalid(
                            "session_terminal.pause_fence",
                            format!(
                                "more than {MAX_SESSION_FLUSH_TERMINALS} Pending submissions require terminalization before one pause acknowledgement"
                            ),
                        ));
                    }
                }
            }
            if !has_more {
                break;
            }
            cursor = last;
        }
        Ok(authorities)
    }

    /// Commit one immutable pre-execution claim.
    ///
    /// A fresh claim returns the only positive terminalization permit. An
    /// existing-identical Pending claim never returns a permit, and a terminal
    /// claim returns the original verified receipt. Reusing an authority with
    /// different identity or payload bytes conflicts atomically.
    ///
    /// # Errors
    /// Foreign ledger bindings, malformed or oversized claims, authority
    /// conflicts, open transactions, corruption, and ledger failures are
    /// structured errors.
    pub fn claim_session_mutation(
        &self,
        claim: &SessionMutationClaim<'_>,
    ) -> Result<SessionMutationClaimResult, LedgerError> {
        if self.in_transaction() {
            return Err(invalid(
                "session_claim.transaction",
                "claim persistence must own its transaction; commit or roll back first",
            ));
        }
        let current_ledger = self.checked_instance_id()?;
        let prepared = validate_claim(*claim, current_ledger)?;
        self.begin()?;
        let write = (|| match self.session_mutation_claim(&claim.authority)? {
            Some(stored) if stored.matches(*claim, prepared) => {
                if let Some(terminal) = self.session_terminal(&claim.authority)? {
                    Ok(SessionMutationClaimResult::Terminal {
                        terminal: Box::new(terminal),
                    })
                } else {
                    Ok(SessionMutationClaimResult::Pending {
                        claim: Box::new(stored),
                    })
                }
            }
            Some(_) => Err(invalid(
                "session_claim.authority",
                format!(
                    "authority {} already stores different claim identity or payload bytes",
                    claim.authority
                ),
            )),
            None => {
                if let (Some(ordinal), Some(owner)) =
                    (claim.causal_ordinal, self.causal_ordinal_owner(*claim)?)
                {
                    return Err(invalid(
                        "session_claim.causal_ordinal",
                        format!(
                            "governor/kind causal ordinal {ordinal} is already owned by {owner}"
                        ),
                    ));
                }
                if claim.kind == PRECLAIM_REQUIRED_SUBMISSION_KIND
                    && let Some(pause_authority) = self.terminalized_pause_successor(*claim)?
                {
                    return Err(invalid(
                        "session_claim.generation",
                        format!(
                            "submission generation {} is already fenced by terminal pause acknowledgement {pause_authority}",
                            claim.generation
                        ),
                    ));
                }
                self.insert_session_claim(*claim, prepared)?;
                Ok(SessionMutationClaimResult::Claimed {
                    permit: SessionClaimPermit {
                        authority: claim.authority,
                        ledger_instance_id: current_ledger,
                        claim_hash: prepared.claim_hash,
                    },
                })
            }
        })();

        match write {
            Ok(result) => {
                if let Err(error) = self.commit() {
                    let _ = self.rollback();
                    return Err(error);
                }
                Ok(result)
            }
            Err(error) => {
                let _ = self.rollback();
                Err(error)
            }
        }
    }

    fn verify_batch_terminals(
        &self,
        batch: &SessionTerminalBatch<'_>,
        prepared: &PreparedBatch,
    ) -> Result<(), LedgerError> {
        let member_rows = self
            .conn
            .query_with_params(
                "SELECT authority FROM session_flush_batch_members \
                 WHERE batch_id = ?1 ORDER BY seq LIMIT 1025",
                &[blob_param(prepared.batch_id.as_bytes())],
            )
            .map_err(|error| sql_err("session flush exact batch-members", &error))?;
        if member_rows.len() != batch.groups.len() {
            return Err(stored_corrupt(
                prepared.batch_id,
                "canonical batch membership count changed before exact replay",
            ));
        }
        for (index, ((group, terminal), member)) in batch
            .groups
            .iter()
            .zip(&prepared.terminals)
            .zip(&member_rows)
            .enumerate()
        {
            let member_authority = ContentHash(row_blob(
                member,
                0,
                prepared.batch_id,
                "batch_member.authority",
            )?);
            if member_authority != group.terminal.claim.authority {
                return Err(stored_corrupt(
                    prepared.batch_id,
                    format!("canonical batch membership differs at sequence {index}"),
                ));
            }
            let stored = self
                .session_terminal(&group.terminal.claim.authority)?
                .ok_or_else(|| {
                    stored_corrupt(
                        group.terminal.claim.authority,
                        "committed flush batch references a missing terminal row",
                    )
                })?;
            if !stored.matches(group.terminal, *terminal) {
                return Err(invalid(
                    "session_terminal.authority",
                    format!(
                        "authority {} already stores different terminal or owned-event bytes",
                        group.terminal.claim.authority
                    ),
                ));
            }
        }
        Ok(())
    }

    /// Atomically exact-insert a canonical terminal/event batch.
    ///
    /// A group with no existing claim may atomically insert claim+terminal
    /// when permit is None, which is the lane for already-completed
    /// non-execution mutations. `submission` is explicitly excluded: it must
    /// already have its pre-execution claim and exact positive permit. Any
    /// other existing Pending claim also requires its exact permit; Pending
    /// plus no permit is Indeterminate and refused.
    /// A new terminal appends each global event, captures its rowid, and
    /// inserts an immutable authority/sequence ownership link in the same
    /// transaction. Existing-identical terminals append nothing. The batch id
    /// is computed internally from the checked ledger and complete ordered
    /// group preimage.
    ///
    /// # Errors
    /// Missing claims, invalid permits, bounds, malformed JSON, authority
    /// conflicts, open transactions, corruption, and ledger failures fail the
    /// complete transaction.
    #[allow(clippy::too_many_lines)] // One explicit all-or-nothing persistence protocol.
    pub fn append_session_terminal_batch(
        &self,
        batch: &SessionTerminalBatch<'_>,
    ) -> Result<SessionTerminalBatchResult, LedgerError> {
        if self.in_transaction() {
            return Err(invalid(
                "session_flush_batch.transaction",
                "batch persistence must own its transaction; commit or roll back first",
            ));
        }
        let prepared = prepare_batch(self, batch)?;
        self.begin()?;
        let write = (|| {
            if let Some(stored_batch) = self.session_flush_batch(&prepared.batch_id)? {
                if !stored_batch.matches(&prepared, batch.groups.len()) {
                    return Err(stored_corrupt(
                        prepared.batch_id,
                        "canonical batch marker disagrees with its derived totals",
                    ));
                }
                self.verify_batch_terminals(batch, &prepared)?;
                return Ok(SessionTerminalBatchResult::Replayed {
                    batch_id: prepared.batch_id,
                });
            }

            let batch_authorities: BTreeSet<_> = batch
                .groups
                .iter()
                .map(|group| group.terminal.claim.authority)
                .collect();
            for group in batch.groups {
                if group.terminal.claim.kind == PAUSE_ACKNOWLEDGEMENT_KIND {
                    let pending = self.pending_submission_predecessors(group.terminal.claim)?;
                    if let Some(authority) = pending
                        .into_iter()
                        .find(|authority| !batch_authorities.contains(authority))
                    {
                        return Err(invalid(
                            "session_terminal.pause_fence",
                            format!(
                                "pause acknowledgement {} cannot terminalize while submission {authority} remains durably Pending in the draining generation",
                                group.terminal.claim.authority
                            ),
                        ));
                    }
                }
            }

            let mut terminals_inserted = 0usize;
            let mut events_appended = 0usize;
            for (group, terminal) in batch.groups.iter().zip(&prepared.terminals) {
                let authority = group.terminal.claim.authority;
                let claim_exists = match self.session_mutation_claim(&authority)? {
                    Some(stored_claim)
                        if stored_claim.matches(group.terminal.claim, terminal.claim) =>
                    {
                        true
                    }
                    Some(_) => {
                        return Err(invalid(
                            "session_claim.authority",
                            format!(
                                "authority {authority} stores different claim identity or payload bytes"
                            ),
                        ));
                    }
                    None => false,
                };
                match self.session_terminal(&authority)? {
                    Some(stored) if stored.matches(group.terminal, *terminal) => {}
                    Some(_) => {
                        return Err(invalid(
                            "session_terminal.authority",
                            format!(
                                "authority {authority} stores different terminal or owned-event bytes"
                            ),
                        ));
                    }
                    None => {
                        if claim_exists && group.terminal.permit.is_none() {
                            return Err(invalid(
                                "session_terminal.claim",
                                format!(
                                    "authority {authority} is durably Pending; without its positive permit the outcome is Indeterminate and work must not be terminalized"
                                ),
                            ));
                        }
                        if !claim_exists {
                            if group.terminal.permit.is_some() {
                                return Err(stored_corrupt(
                                    authority,
                                    "positive claim permit exists but its durable claim is missing",
                                ));
                            }
                            if group.terminal.claim.kind == PRECLAIM_REQUIRED_SUBMISSION_KIND {
                                return Err(invalid(
                                    "session_terminal.claim",
                                    format!(
                                        "submission authority {authority} requires a committed pre-execution claim and its positive terminalization permit"
                                    ),
                                ));
                            }
                            self.insert_session_claim(group.terminal.claim, terminal.claim)?;
                        }
                        self.insert_session_terminal(group.terminal, *terminal)?;
                        terminals_inserted += 1;
                        for (seq, event) in group.events.iter().enumerate() {
                            let event_id = self.append_event(event)?;
                            self.insert_session_terminal_event_link(authority, seq, event_id)?;
                            events_appended += 1;
                        }
                    }
                }
            }
            for group in batch.groups {
                let authority = group.terminal.claim.authority;
                let memberships = self.terminal_batch_memberships(&authority)?;
                if memberships.len() >= MAX_SESSION_TERMINAL_BATCH_MEMBERSHIPS {
                    return Err(invalid(
                        "session_terminal.batch_memberships",
                        format!(
                            "terminal {authority} already has the {MAX_SESSION_TERMINAL_BATCH_MEMBERSHIPS}-batch witness limit"
                        ),
                    ));
                }
            }
            self.insert_session_flush_batch(&prepared, batch.groups.len())?;
            for (seq, group) in batch.groups.iter().enumerate() {
                self.insert_session_flush_batch_member(
                    prepared.batch_id,
                    seq,
                    group.terminal.claim.authority,
                )?;
            }
            Ok(SessionTerminalBatchResult::Committed {
                batch_id: prepared.batch_id,
                terminals_inserted,
                events_appended,
            })
        })();

        match write {
            Ok(result) => {
                if let Err(error) = self.commit() {
                    let _ = self.rollback();
                    return Err(error);
                }
                Ok(result)
            }
            Err(error) => {
                let _ = self.rollback();
                Err(error)
            }
        }
    }
}

fn hex_bytes(bytes: &[u8]) -> String {
    use core::fmt::Write as _;

    let mut encoded = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        let _ = write!(encoded, "{byte:02x}");
    }
    encoded
}

#[cfg(test)]
mod tests {
    use super::*;

    fn authority(seed: u64) -> ContentHash {
        hash_bytes(&seed.to_le_bytes())
    }

    fn drained_report(run: fs_exec::RunId) -> fs_exec::DrainFinalizeReport {
        let gate = fs_exec::CancelGate::new_clock_free();
        let tracker = fs_exec::DrainTracker::new(run, &gate);
        let worker = tracker.register_worker().expect("fixture worker");
        gate.request();
        drop(worker);
        tracker.finalize().expect("fixture run drained")
    }

    #[test]
    fn solver_checkpoint_receipt_is_idempotent_transport_checked_and_ledger_bound() {
        let ledger = Ledger::open(":memory:").expect("checkpoint ledger");
        let run = fs_exec::RunId(0xA771);
        let report = drained_report(run);
        let snapshot =
            fs_exec::solver::envelope::seal(0x4653_434b_5054, 1, run.0, b"bounded-solver-state");
        let artifact = ledger
            .put_artifact(SOLVER_STATE_ARTIFACT_KIND, &snapshot, None)
            .expect("solver-state artifact");
        let pause_authority = authority(0x51);
        let claim = SolverCheckpointClaim {
            session: 71,
            pause_authority,
            gate_generation: 4,
            solver_state_artifact: artifact.hash,
            drain_report: &report,
        };
        let first = ledger
            .attest_solver_checkpoint(claim)
            .expect("first attestation");
        let after_response_loss = ledger
            .attest_solver_checkpoint(claim)
            .expect("exact retry after response loss");
        assert_eq!(after_response_loss, first);
        assert_eq!(
            ledger
                .solver_checkpoint_receipt(pause_authority)
                .expect("receipt lookup"),
            Some(first)
        );

        let transport = first.to_bytes();
        let candidate = SolverCheckpointReceipt::from_bytes(&transport).expect("typed transport");
        ledger
            .verify_solver_checkpoint_receipt(&candidate)
            .expect("stored transport candidate verifies");
        let mut forged = transport;
        forged[20] ^= 1;
        assert!(SolverCheckpointReceipt::from_bytes(&forged).is_err());

        let foreign = Ledger::open(":memory:").expect("foreign ledger");
        assert!(foreign.verify_solver_checkpoint_receipt(&first).is_err());

        let changed_snapshot =
            fs_exec::solver::envelope::seal(0x4653_434b_5054, 1, run.0, b"changed-solver-state");
        let changed_artifact = ledger
            .put_artifact(SOLVER_STATE_ARTIFACT_KIND, &changed_snapshot, None)
            .expect("changed solver-state artifact");
        let conflicting = SolverCheckpointClaim {
            solver_state_artifact: changed_artifact.hash,
            ..claim
        };
        assert!(ledger.attest_solver_checkpoint(conflicting).is_err());
        assert_eq!(
            ledger
                .solver_checkpoint_receipt(pause_authority)
                .expect("original receipt survives conflict"),
            Some(first)
        );
    }

    #[test]
    fn solver_checkpoint_attestation_refuses_wrong_kind_and_run() {
        let ledger = Ledger::open(":memory:").expect("checkpoint refusal ledger");
        let run = fs_exec::RunId(0xA772);
        let report = drained_report(run);
        let wrong_run = fs_exec::solver::envelope::seal(7, 1, run.0 + 1, b"state");
        let wrong_run_artifact = ledger
            .put_artifact(SOLVER_STATE_ARTIFACT_KIND, &wrong_run, None)
            .expect("wrong-run artifact");
        assert!(
            ledger
                .attest_solver_checkpoint(SolverCheckpointClaim {
                    session: 72,
                    pause_authority: authority(0x52),
                    gate_generation: 0,
                    solver_state_artifact: wrong_run_artifact.hash,
                    drain_report: &report,
                })
                .is_err()
        );

        let valid_envelope = fs_exec::solver::envelope::seal(7, 1, run.0, b"state");
        let wrong_kind = ledger
            .put_artifact("generic-artifact", &valid_envelope, None)
            .expect("wrong-kind artifact");
        assert!(
            ledger
                .attest_solver_checkpoint(SolverCheckpointClaim {
                    session: 72,
                    pause_authority: authority(0x53),
                    gate_generation: 0,
                    solver_state_artifact: wrong_kind.hash,
                    drain_report: &report,
                })
                .is_err()
        );
    }

    fn fixture_claim<'a>(
        ledger: &Ledger,
        authority: ContentHash,
        payload: &'a [u8],
    ) -> SessionMutationClaim<'a> {
        SessionMutationClaim {
            authority,
            ledger_instance_id: ledger.instance_id(),
            governor_hash: hash_bytes(b"registry-corruption-governor"),
            session_open_hash: hash_bytes(b"registry-corruption-open"),
            kind: "test-atomic",
            session: 71,
            ledger_scope: "registry-corruption",
            generation: 2,
            causal_ordinal: None,
            payload,
        }
    }

    fn commit_fixture(ledger: &Ledger, authority: ContentHash) {
        let claim = fixture_claim(ledger, authority, b"payload");
        let session = 71_u64.to_be_bytes();
        let event = EventRow {
            session: Some(&session),
            t: 9,
            kind: "session.idempotent-execution",
            payload: Some(r#"{"schema":"registry-corruption-v1"}"#),
        };
        let events = [event];
        let group = SessionTerminalGroup {
            terminal: SessionTerminalRow {
                claim,
                permit: None,
                receipt: b"receipt",
            },
            events: &events,
        };
        let groups = [group];
        ledger
            .append_session_terminal_batch(&SessionTerminalBatch { groups: &groups })
            .expect("fixture terminal batch");
    }

    /// Insert hash-valid registry rows inside a caller-owned test transaction.
    /// This bypasses only public admission so exact read caps can be exercised
    /// without thousands of nested transactions.
    fn insert_canonical_terminal_batch_fixture(
        ledger: &Ledger,
        claims: &[SessionMutationClaim<'_>],
    ) {
        let groups: Vec<_> = claims
            .iter()
            .copied()
            .map(|claim| SessionTerminalGroup {
                terminal: SessionTerminalRow {
                    claim,
                    permit: None,
                    receipt: b"terminal",
                },
                events: &[],
            })
            .collect();
        let batch = SessionTerminalBatch { groups: &groups };
        let prepared = prepare_batch(ledger, &batch).expect("canonical fixture batch");
        for (group, terminal) in groups.iter().zip(&prepared.terminals) {
            ledger
                .insert_session_claim(group.terminal.claim, terminal.claim)
                .expect("canonical fixture claim");
            ledger
                .insert_session_terminal(group.terminal, *terminal)
                .expect("canonical fixture terminal");
        }
        ledger
            .insert_session_flush_batch(&prepared, groups.len())
            .expect("canonical fixture batch marker");
        for (seq, group) in groups.iter().enumerate() {
            ledger
                .insert_session_flush_batch_member(
                    prepared.batch_id,
                    seq,
                    group.terminal.claim.authority,
                )
                .expect("canonical fixture batch member");
        }
    }

    /// Insert the production-shaped common submission recovery row: one
    /// submission claim, one canonical owned event, and one singleton batch.
    fn insert_canonical_submission_terminal_fixture(
        ledger: &Ledger,
        claim: SessionMutationClaim<'_>,
    ) {
        assert_eq!(claim.kind, PRECLAIM_REQUIRED_SUBMISSION_KIND);
        let event_ordinal = claim
            .causal_ordinal
            .expect("submission fixture carries an admission ordinal");
        let session = claim.session.to_be_bytes();
        let event = EventRow {
            session: Some(&session),
            t: i64::try_from(event_ordinal).expect("bounded submission ordinal"),
            kind: "session.idempotent-execution",
            payload: Some(r#"{"schema":"fs-session-idempotency-v5","result":"done"}"#),
        };
        let events = [event];
        let groups = [SessionTerminalGroup {
            terminal: SessionTerminalRow {
                claim,
                permit: None,
                receipt: b"submission-terminal-v1",
            },
            events: &events,
        }];
        let batch = SessionTerminalBatch { groups: &groups };
        let prepared = prepare_batch(ledger, &batch).expect("canonical submission fixture batch");
        ledger
            .insert_session_claim(claim, prepared.terminals[0].claim)
            .expect("canonical submission fixture claim");
        ledger
            .insert_session_terminal(groups[0].terminal, prepared.terminals[0])
            .expect("canonical submission fixture terminal");
        let event_id = ledger
            .append_event(&event)
            .expect("canonical submission fixture event");
        ledger
            .insert_session_terminal_event_link(claim.authority, 0, event_id)
            .expect("canonical submission fixture event link");
        ledger
            .insert_session_flush_batch(&prepared, 1)
            .expect("canonical submission fixture batch marker");
        ledger
            .insert_session_flush_batch_member(prepared.batch_id, 0, claim.authority)
            .expect("canonical submission fixture batch member");
    }

    /// Add a new two-member witness for one existing anchor terminal plus one
    /// fresh companion, again inside the caller-owned test transaction.
    fn insert_canonical_companion_witness_fixture(
        ledger: &Ledger,
        anchor: SessionMutationClaim<'_>,
        companion: SessionMutationClaim<'_>,
    ) {
        let groups = [
            SessionTerminalGroup {
                terminal: SessionTerminalRow {
                    claim: anchor,
                    permit: None,
                    receipt: b"terminal",
                },
                events: &[],
            },
            SessionTerminalGroup {
                terminal: SessionTerminalRow {
                    claim: companion,
                    permit: None,
                    receipt: b"terminal",
                },
                events: &[],
            },
        ];
        let batch = SessionTerminalBatch { groups: &groups };
        let prepared = prepare_batch(ledger, &batch).expect("canonical witness batch");
        ledger
            .insert_session_claim(companion, prepared.terminals[1].claim)
            .expect("canonical companion claim");
        ledger
            .insert_session_terminal(groups[1].terminal, prepared.terminals[1])
            .expect("canonical companion terminal");
        ledger
            .insert_session_flush_batch(&prepared, groups.len())
            .expect("canonical witness batch marker");
        for (seq, group) in groups.iter().enumerate() {
            ledger
                .insert_session_flush_batch_member(
                    prepared.batch_id,
                    seq,
                    group.terminal.claim.authority,
                )
                .expect("canonical witness batch member");
        }
    }

    #[test]
    fn immutable_registry_and_owned_event_triggers_refuse_all_updates() {
        let ledger = Ledger::open(":memory:").expect("fixture ledger");
        let authority = authority(1);
        commit_fixture(&ledger, authority);
        for statement in [
            "UPDATE session_claims SET payload = X'00'",
            "UPDATE session_claim_discovery SET kind = 'changed'",
            "UPDATE session_terminals SET receipt = X'00'",
            "UPDATE session_terminal_events SET seq = 1",
            "UPDATE session_flush_batches SET terminal_count = 2",
            "UPDATE session_flush_batch_members SET seq = 1",
            "UPDATE events SET payload = '{}'",
        ] {
            assert!(
                ledger.conn.execute(statement).is_err(),
                "immutable write unexpectedly succeeded: {statement}"
            );
        }
        assert!(ledger.session_terminal(&authority).unwrap().is_some());
    }

    #[test]
    fn dual_discovery_witness_keeps_filtered_corruption_fail_closed() {
        let claim_ledger = Ledger::open(":memory:").expect("claim-corruption ledger");
        let claim_authority = authority(11);
        let claim = fixture_claim(&claim_ledger, claim_authority, b"pending");
        assert!(matches!(
            claim_ledger
                .claim_session_mutation(&claim)
                .expect("pending claim fixture"),
            SessionMutationClaimResult::Claimed { .. }
        ));
        claim_ledger
            .conn
            .execute("DROP TRIGGER trg_session_claims_immutable_update")
            .expect("test-only claim trigger bypass");
        claim_ledger
            .conn
            .prepare("UPDATE session_claims SET governor_hash = ?1 WHERE authority = ?2")
            .expect("prepare hidden-claim corruption")
            .execute_with_params(&[
                blob_param(hash_bytes(b"foreign-governor").as_bytes()),
                blob_param(claim_authority.as_bytes()),
            ])
            .expect("inject hidden-claim corruption");
        assert!(matches!(
            claim_ledger.session_mutation_claim_count(claim.governor_hash),
            Err(LedgerError::Corrupt { .. })
        ));
        assert!(matches!(
            claim_ledger.pending_session_mutation(
                claim.governor_hash,
                claim.session_open_hash,
                claim.kind,
                claim.session,
                claim.ledger_scope,
                claim.generation,
            ),
            Err(LedgerError::Corrupt { .. })
        ));

        let witness_ledger = Ledger::open(":memory:").expect("witness-corruption ledger");
        let witness_authority = authority(12);
        let witness_claim = fixture_claim(&witness_ledger, witness_authority, b"pending");
        assert!(matches!(
            witness_ledger
                .claim_session_mutation(&witness_claim)
                .expect("witness claim fixture"),
            SessionMutationClaimResult::Claimed { .. }
        ));
        witness_ledger
            .conn
            .execute("DROP TRIGGER trg_session_claim_discovery_immutable_update")
            .expect("test-only discovery trigger bypass");
        witness_ledger
            .conn
            .prepare("UPDATE session_claim_discovery SET ledger_scope = ?1 WHERE authority = ?2")
            .expect("prepare hidden-witness corruption")
            .execute_with_params(&[
                text_param("different-scope"),
                blob_param(witness_authority.as_bytes()),
            ])
            .expect("inject hidden-witness corruption");
        assert!(matches!(
            witness_ledger.pending_session_mutation(
                witness_claim.governor_hash,
                witness_claim.session_open_hash,
                witness_claim.kind,
                witness_claim.session,
                witness_claim.ledger_scope,
                witness_claim.generation,
            ),
            Err(LedgerError::Corrupt { .. })
        ));

        let missing_ledger = Ledger::open(":memory:").expect("missing-witness ledger");
        let missing_authority = authority(13);
        let missing_claim = fixture_claim(&missing_ledger, missing_authority, b"pending");
        assert!(matches!(
            missing_ledger
                .claim_session_mutation(&missing_claim)
                .expect("missing-witness claim fixture"),
            SessionMutationClaimResult::Claimed { .. }
        ));
        missing_ledger
            .conn
            .execute("DROP TRIGGER trg_session_claim_discovery_immutable_delete")
            .expect("test-only discovery delete trigger bypass");
        missing_ledger
            .conn
            .prepare("DELETE FROM session_claim_discovery WHERE authority = ?1")
            .expect("prepare missing-witness corruption")
            .execute_with_params(&[blob_param(missing_authority.as_bytes())])
            .expect("inject missing-witness corruption");
        assert!(matches!(
            missing_ledger.session_mutation_claim(&missing_authority),
            Err(LedgerError::Corrupt { .. })
        ));
    }

    #[test]
    fn future_claim_schema_and_receipt_hash_tampering_fail_closed() {
        let future = Ledger::open(":memory:").expect("future-schema ledger");
        let future_authority = authority(2);
        commit_fixture(&future, future_authority);
        future
            .conn
            .execute("ALTER TABLE session_claims RENAME TO session_claims_v1_fixture")
            .expect("move checked table behind a test-only name");
        future
            .conn
            .execute("CREATE TABLE session_claims AS SELECT * FROM session_claims_v1_fixture")
            .expect("copy exact rows into a test-only unconstrained shadow");
        future
            .conn
            .prepare("UPDATE session_claims SET registry_schema_version = 2 WHERE authority = ?1")
            .unwrap()
            .execute_with_params(&[blob_param(future_authority.as_bytes())])
            .expect("inject future registry envelope");
        assert!(matches!(
            future.session_mutation_claim(&future_authority),
            Err(LedgerError::Corrupt { .. })
        ));

        let tampered = Ledger::open(":memory:").expect("receipt-tamper ledger");
        let tampered_authority = authority(3);
        commit_fixture(&tampered, tampered_authority);
        tampered
            .conn
            .execute("DROP TRIGGER trg_session_terminals_immutable_update")
            .expect("test-only trigger bypass");
        tampered
            .conn
            .prepare("UPDATE session_terminals SET receipt = X'00' WHERE authority = ?1")
            .unwrap()
            .execute_with_params(&[blob_param(tampered_authority.as_bytes())])
            .expect("inject receipt/hash mismatch");
        assert!(matches!(
            tampered.session_terminal(&tampered_authority),
            Err(LedgerError::Corrupt { .. })
        ));
    }

    #[test]
    fn missing_event_link_or_batch_member_is_detected_before_replay() {
        let ledger = Ledger::open(":memory:").expect("fixture ledger");
        let event_authority = authority(4);
        commit_fixture(&ledger, event_authority);
        ledger
            .conn
            .execute("DROP TRIGGER trg_session_terminal_events_immutable_delete")
            .expect("test-only event-link trigger bypass");
        ledger
            .conn
            .prepare("DELETE FROM session_terminal_events WHERE authority = ?1")
            .unwrap()
            .execute_with_params(&[blob_param(event_authority.as_bytes())])
            .expect("inject partial owned-event group");
        assert!(matches!(
            ledger.session_terminal(&event_authority),
            Err(LedgerError::Corrupt { .. })
        ));

        let batch_ledger = Ledger::open(":memory:").expect("batch-member ledger");
        let batch_authority = authority(5);
        commit_fixture(&batch_ledger, batch_authority);
        batch_ledger
            .conn
            .execute("DROP TRIGGER trg_session_flush_batch_members_immutable_delete")
            .expect("test-only batch-member trigger bypass");
        batch_ledger
            .conn
            .execute("DELETE FROM session_flush_batch_members")
            .expect("inject partial batch membership");
        let claim = fixture_claim(&batch_ledger, batch_authority, b"payload");
        let session = 71_u64.to_be_bytes();
        let event = EventRow {
            session: Some(&session),
            t: 9,
            kind: "session.idempotent-execution",
            payload: Some(r#"{"schema":"registry-corruption-v1"}"#),
        };
        let events = [event];
        let group = SessionTerminalGroup {
            terminal: SessionTerminalRow {
                claim,
                permit: None,
                receipt: b"receipt",
            },
            events: &events,
        };
        let groups = [group];
        assert!(matches!(
            batch_ledger.append_session_terminal_batch(&SessionTerminalBatch { groups: &groups }),
            Err(LedgerError::Corrupt { .. })
        ));

        let marker_ledger = Ledger::open(":memory:").expect("batch-marker ledger");
        let marker_authority = authority(6);
        commit_fixture(&marker_ledger, marker_authority);
        marker_ledger
            .conn
            .execute("DROP TRIGGER trg_session_flush_batches_immutable_update")
            .expect("test-only batch-marker trigger bypass");
        marker_ledger
            .conn
            .execute("UPDATE session_flush_batches SET event_count = 0")
            .expect("inject inconsistent batch-marker totals");
        assert!(matches!(
            marker_ledger.session_terminal(&marker_authority),
            Err(LedgerError::Corrupt { .. })
        ));
    }

    #[test]
    fn simple_generation_probe_handles_pending_and_one_event_terminal() {
        let pending_ledger = Ledger::open(":memory:").expect("pending-probe ledger");
        let pending_authority = authority(7);
        let pending_claim = fixture_claim(&pending_ledger, pending_authority, b"pending");
        assert!(matches!(
            pending_ledger
                .claim_session_mutation(&pending_claim)
                .expect("pending fixture claim"),
            SessionMutationClaimResult::Claimed { .. }
        ));
        let reads_before = pending_ledger.read_queries();
        let stored = pending_ledger
            .pending_session_mutation(
                pending_claim.governor_hash,
                pending_claim.session_open_hash,
                pending_claim.kind,
                pending_claim.session,
                pending_claim.ledger_scope,
                pending_claim.generation,
            )
            .expect("simple Pending probe")
            .expect("Pending claim");
        assert_eq!(stored.authority, pending_authority);
        assert_eq!(stored.payload.as_slice(), pending_claim.payload);
        assert_eq!(
            pending_ledger.read_queries() - reads_before,
            2,
            "simple Pending recovery must not recurse through typed readers"
        );

        let event_ledger = Ledger::open(":memory:").expect("event-probe ledger");
        let event_authority = authority(8);
        let event_claim = SessionMutationClaim {
            kind: PRECLAIM_REQUIRED_SUBMISSION_KIND,
            causal_ordinal: Some(1),
            ..fixture_claim(&event_ledger, event_authority, b"submission")
        };
        event_ledger.begin().expect("one-event fixture transaction");
        insert_canonical_submission_terminal_fixture(&event_ledger, event_claim);
        event_ledger.commit().expect("one-event fixture commit");
        let probe = event_ledger
            .simple_generation_probe(event_authority)
            .expect("one-event probe")
            .expect("canonical one-event fast path");
        assert_eq!(probe.claim.authority, event_authority);
        assert!(probe.terminalized);
        let reads_before = event_ledger.read_queries();
        assert_eq!(
            event_ledger
                .pending_session_mutation(
                    event_claim.governor_hash,
                    event_claim.session_open_hash,
                    event_claim.kind,
                    event_claim.session,
                    event_claim.ledger_scope,
                    event_claim.generation,
                )
                .expect("one-event recovery probe"),
            None
        );
        assert_eq!(
            event_ledger.read_queries() - reads_before,
            2,
            "one-event recovery must not recurse through typed readers"
        );
    }

    #[test]
    fn simple_generation_probe_defers_malformed_owned_event() {
        let ledger = Ledger::open(":memory:").expect("malformed-event ledger");
        let authority = authority(10);
        let claim = SessionMutationClaim {
            kind: PRECLAIM_REQUIRED_SUBMISSION_KIND,
            causal_ordinal: Some(1),
            ..fixture_claim(&ledger, authority, b"submission")
        };
        ledger.begin().expect("malformed-event fixture transaction");
        insert_canonical_submission_terminal_fixture(&ledger, claim);
        ledger.commit().expect("malformed-event fixture commit");
        ledger
            .conn
            .execute("DROP TRIGGER trg_owned_session_events_immutable_update")
            .expect("test-only owned-event trigger bypass");
        ledger
            .conn
            .execute("UPDATE events SET kind = 'bad kind'")
            .expect("inject non-canonical event kind");

        assert!(
            ledger
                .simple_generation_probe(authority)
                .expect("malformed event probe")
                .is_none(),
            "malformed events must defer to the complete verifier"
        );
        assert!(matches!(
            ledger.pending_session_mutation(
                claim.governor_hash,
                claim.session_open_hash,
                claim.kind,
                claim.session,
                claim.ledger_scope,
                claim.generation,
            ),
            Err(LedgerError::Corrupt { .. })
        ));
    }

    #[test]
    fn simple_generation_probe_masks_oversized_terminal_scalars() {
        let ledger = Ledger::open(":memory:").expect("scalar-mask ledger");
        let authority = authority(9);
        let claim = fixture_claim(&ledger, authority, b"scalar-mask");
        ledger.begin().expect("scalar-mask fixture transaction");
        insert_canonical_terminal_batch_fixture(&ledger, &[claim]);
        ledger.commit().expect("scalar-mask fixture commit");

        ledger
            .conn
            .execute("ALTER TABLE session_terminals RENAME TO session_terminals_scalar_fixture")
            .expect("move checked terminal table behind a test-only name");
        ledger
            .conn
            .execute(
                "CREATE TABLE session_terminals AS \
                 SELECT * FROM session_terminals_scalar_fixture",
            )
            .expect("copy terminal rows into a test-only unconstrained shadow");
        let oversized_scalar = vec![0_u8; MAX_SESSION_FLUSH_ENCODED_BYTES + 1];
        ledger
            .conn
            .prepare("UPDATE session_terminals SET event_count = ?1 WHERE authority = ?2")
            .expect("prepare oversized scalar injection")
            .execute_with_params(&[
                blob_param(&oversized_scalar),
                blob_param(authority.as_bytes()),
            ])
            .expect("inject oversized scalar storage");

        assert!(
            ledger
                .simple_generation_probe(authority)
                .expect("bounded simple probe")
                .is_none(),
            "malformed terminal scalars must defer to the deep verifier"
        );
        assert!(matches!(
            ledger.pending_session_mutation(
                claim.governor_hash,
                claim.session_open_hash,
                claim.kind,
                claim.session,
                claim.ledger_scope,
                claim.generation,
            ),
            Err(LedgerError::Corrupt { .. })
        ));
    }

    #[test]
    fn recovery_probe_accepts_exact_claim_cap_and_rejects_limit_plus_one() {
        let ledger = Ledger::open(":memory:").expect("recovery-cap ledger");
        ledger.begin().expect("bulk fixture transaction");
        for index in 0..MAX_SESSION_RECOVERY_PROBE_CLAIMS {
            let claim = SessionMutationClaim {
                kind: PRECLAIM_REQUIRED_SUBMISSION_KIND,
                causal_ordinal: Some(u64::try_from(index + 1).expect("bounded ordinal")),
                ..fixture_claim(
                    &ledger,
                    authority(100_000 + u64::try_from(index).expect("bounded index")),
                    b"",
                )
            };
            insert_canonical_submission_terminal_fixture(&ledger, claim);
        }
        ledger.commit().expect("exact-cap fixture commit");

        let envelope = SessionMutationClaim {
            kind: PRECLAIM_REQUIRED_SUBMISSION_KIND,
            causal_ordinal: Some(1),
            ..fixture_claim(&ledger, authority(100_000), b"")
        };
        let reads_before = ledger.read_queries();
        assert_eq!(
            ledger
                .pending_session_mutation(
                    envelope.governor_hash,
                    envelope.session_open_hash,
                    envelope.kind,
                    envelope.session,
                    envelope.ledger_scope,
                    envelope.generation,
                )
                .expect("exact recovery-probe cap"),
            None
        );
        assert_eq!(
            ledger.read_queries() - reads_before,
            2,
            "the exact-cap recovery probe must not recurse through per-claim typed readers"
        );

        ledger.begin().expect("limit-plus-one transaction");
        let overflow = SessionMutationClaim {
            kind: PRECLAIM_REQUIRED_SUBMISSION_KIND,
            causal_ordinal: Some(
                u64::try_from(MAX_SESSION_RECOVERY_PROBE_CLAIMS + 1).expect("bounded ordinal"),
            ),
            ..fixture_claim(
                &ledger,
                authority(
                    100_000
                        + u64::try_from(MAX_SESSION_RECOVERY_PROBE_CLAIMS).expect("bounded cap"),
                ),
                b"",
            )
        };
        insert_canonical_submission_terminal_fixture(&ledger, overflow);
        ledger.commit().expect("limit-plus-one fixture commit");
        assert_eq!(
            ledger.pending_session_mutation(
                envelope.governor_hash,
                envelope.session_open_hash,
                envelope.kind,
                envelope.session,
                envelope.ledger_scope,
                envelope.generation,
            ),
            Err(LedgerError::Invalid {
                field: "session_claim.recovery_probe".to_string(),
                problem: format!(
                    "generation contains more than the {MAX_SESSION_RECOVERY_PROBE_CLAIMS}-claim recovery-probe limit"
                ),
            })
        );
    }

    #[test]
    fn pause_fence_accepts_exact_submission_cap_and_rejects_limit_plus_one() {
        let ledger = Ledger::open(":memory:").expect("pause-fence-cap ledger");
        ledger.begin().expect("bulk fixture transaction");
        for index in 0..MAX_SESSION_PAUSE_FENCE_SUBMISSIONS {
            let claim = SessionMutationClaim {
                kind: PRECLAIM_REQUIRED_SUBMISSION_KIND,
                causal_ordinal: Some(u64::try_from(index + 1).expect("bounded ordinal")),
                ..fixture_claim(
                    &ledger,
                    authority(200_000 + u64::try_from(index).expect("bounded index")),
                    b"",
                )
            };
            insert_canonical_submission_terminal_fixture(&ledger, claim);
        }
        ledger.commit().expect("exact-cap fixture commit");

        let pause = SessionMutationClaim {
            authority: authority(210_000),
            kind: PAUSE_ACKNOWLEDGEMENT_KIND,
            generation: 3,
            causal_ordinal: Some(1),
            payload: b"pause",
            ..fixture_claim(&ledger, authority(210_000), b"pause")
        };
        let reads_before = ledger.read_queries();
        assert!(
            ledger
                .pending_submission_predecessors(pause)
                .expect("exact pause-fence cap")
                .is_empty()
        );
        assert_eq!(
            ledger.read_queries() - reads_before,
            0,
            "the exact-cap pause fence must not recurse through per-claim typed readers"
        );

        ledger.begin().expect("limit-plus-one transaction");
        let overflow_index = MAX_SESSION_PAUSE_FENCE_SUBMISSIONS;
        let overflow = SessionMutationClaim {
            kind: PRECLAIM_REQUIRED_SUBMISSION_KIND,
            causal_ordinal: Some(u64::try_from(overflow_index + 1).expect("bounded ordinal")),
            ..fixture_claim(
                &ledger,
                authority(200_000 + u64::try_from(overflow_index).expect("bounded index")),
                b"",
            )
        };
        insert_canonical_submission_terminal_fixture(&ledger, overflow);
        ledger.commit().expect("limit-plus-one fixture commit");
        assert_eq!(
            ledger.pending_submission_predecessors(pause),
            Err(LedgerError::Invalid {
                field: "session_terminal.pause_fence".to_string(),
                problem: format!(
                    "draining generation contains more than the {MAX_SESSION_PAUSE_FENCE_SUBMISSIONS}-submission fence limit"
                ),
            })
        );
    }

    #[test]
    fn terminal_read_accepts_exact_membership_cap_and_rejects_limit_plus_one() {
        let ledger = Ledger::open(":memory:").expect("membership-cap ledger");
        let anchor = fixture_claim(&ledger, authority(300_000), b"anchor");
        ledger.begin().expect("bulk fixture transaction");
        insert_canonical_terminal_batch_fixture(&ledger, &[anchor]);
        for index in 1..MAX_SESSION_TERMINAL_BATCH_MEMBERSHIPS {
            let companion = fixture_claim(
                &ledger,
                authority(300_000 + u64::try_from(index).expect("bounded index")),
                b"companion",
            );
            insert_canonical_companion_witness_fixture(&ledger, anchor, companion);
        }
        ledger.commit().expect("exact-cap fixture commit");

        assert_eq!(
            ledger
                .terminal_batch_memberships(&anchor.authority)
                .expect("exact membership cap")
                .len(),
            MAX_SESSION_TERMINAL_BATCH_MEMBERSHIPS
        );
        assert!(
            ledger
                .session_terminal(&anchor.authority)
                .expect("exact-cap terminal read")
                .is_some()
        );

        ledger.begin().expect("limit-plus-one transaction");
        let overflow_companion = fixture_claim(
            &ledger,
            authority(
                300_000
                    + u64::try_from(MAX_SESSION_TERMINAL_BATCH_MEMBERSHIPS).expect("bounded cap"),
            ),
            b"companion",
        );
        insert_canonical_companion_witness_fixture(&ledger, anchor, overflow_companion);
        ledger.commit().expect("limit-plus-one fixture commit");
        assert_eq!(
            ledger.session_terminal(&anchor.authority),
            Err(LedgerError::Corrupt {
                hash_hex: anchor.authority.to_hex(),
                detail: format!(
                    "session mutation registry: terminal authority exceeds the {MAX_SESSION_TERMINAL_BATCH_MEMBERSHIPS}-batch witness limit"
                ),
            })
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)] // Two independent corrupt-ledger fixtures cover both reciprocal fences.
    fn generation_fences_reject_corrupt_terminal_presence() {
        let predecessor_ledger = Ledger::open(":memory:").expect("predecessor ledger");
        let submission_authority = authority(20);
        let submission = SessionMutationClaim {
            kind: PRECLAIM_REQUIRED_SUBMISSION_KIND,
            causal_ordinal: Some(1),
            ..fixture_claim(&predecessor_ledger, submission_authority, b"submission")
        };
        let permit = match predecessor_ledger
            .claim_session_mutation(&submission)
            .expect("submission claim")
        {
            SessionMutationClaimResult::Claimed { permit } => permit,
            other => panic!("fresh submission returned {other:?}"),
        };
        let submission_groups = [SessionTerminalGroup {
            terminal: SessionTerminalRow {
                claim: submission,
                permit: Some(permit),
                receipt: b"submission-terminal",
            },
            events: &[],
        }];
        predecessor_ledger
            .append_session_terminal_batch(&SessionTerminalBatch {
                groups: &submission_groups,
            })
            .expect("submission terminal");
        predecessor_ledger
            .conn
            .execute("DROP TRIGGER trg_session_flush_batch_members_immutable_delete")
            .expect("test-only batch-member trigger bypass");
        predecessor_ledger
            .conn
            .execute("DELETE FROM session_flush_batch_members")
            .expect("inject corrupt predecessor terminal witness");
        assert!(matches!(
            predecessor_ledger.pending_session_mutation(
                submission.governor_hash,
                submission.session_open_hash,
                submission.kind,
                submission.session,
                submission.ledger_scope,
                submission.generation,
            ),
            Err(LedgerError::Corrupt { .. })
        ));
        let pause = SessionMutationClaim {
            authority: authority(21),
            kind: PAUSE_ACKNOWLEDGEMENT_KIND,
            generation: 3,
            causal_ordinal: Some(2),
            payload: b"pause",
            ..fixture_claim(&predecessor_ledger, authority(21), b"pause")
        };
        let pause_groups = [SessionTerminalGroup {
            terminal: SessionTerminalRow {
                claim: pause,
                permit: None,
                receipt: b"pause-terminal",
            },
            events: &[],
        }];
        assert!(matches!(
            predecessor_ledger.append_session_terminal_batch(&SessionTerminalBatch {
                groups: &pause_groups,
            }),
            Err(LedgerError::Corrupt { .. })
        ));

        let successor_ledger = Ledger::open(":memory:").expect("successor ledger");
        let pause = SessionMutationClaim {
            authority: authority(22),
            kind: PAUSE_ACKNOWLEDGEMENT_KIND,
            generation: 3,
            causal_ordinal: Some(1),
            payload: b"pause",
            ..fixture_claim(&successor_ledger, authority(22), b"pause")
        };
        let pause_groups = [SessionTerminalGroup {
            terminal: SessionTerminalRow {
                claim: pause,
                permit: None,
                receipt: b"pause-terminal",
            },
            events: &[],
        }];
        successor_ledger
            .append_session_terminal_batch(&SessionTerminalBatch {
                groups: &pause_groups,
            })
            .expect("pause terminal");
        successor_ledger
            .conn
            .execute("DROP TRIGGER trg_session_flush_batch_members_immutable_delete")
            .expect("test-only batch-member trigger bypass");
        successor_ledger
            .conn
            .execute("DELETE FROM session_flush_batch_members")
            .expect("inject corrupt successor terminal witness");
        let late_submission = SessionMutationClaim {
            kind: PRECLAIM_REQUIRED_SUBMISSION_KIND,
            causal_ordinal: Some(2),
            ..fixture_claim(&successor_ledger, authority(23), b"late-submission")
        };
        assert!(matches!(
            successor_ledger.claim_session_mutation(&late_submission),
            Err(LedgerError::Corrupt { .. })
        ));
    }
}

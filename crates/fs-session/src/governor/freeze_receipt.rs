//! Predecessor-bound snapshot freeze receipts (bead frankensim-sj31i.52.5.4.1,
//! session lane).
//!
//! The fs-exec side of the bead (crates/fs-exec/src/freeze.rs) produces a
//! committed [`fs_exec::freeze::SnapshotFreezeReceipt`] whose identity binds
//! the owner/session/solver-instance generation, the drained run's executor
//! report, the pause labels, the payload commitment, and the sealed envelope
//! content id. This module is the session-owner half: it validates that such
//! a receipt really belongs to THIS governor and THIS paused generation, then
//! publishes one immutable terminal row plus one event through the same
//! atomic batch machinery every other governor terminal uses, and arms the
//! resume-activation gate so the paused generation cannot run again until the
//! receipt exists.
//!
//! # No-claim boundary
//!
//! A recorded receipt proves the binding checks passed and the batch
//! committed. It does not authenticate the fs-exec registry beyond its
//! identity fields, and it does not validate solver physics.

use crate::SessionError;
use fs_blake3::ContentHash;

pub(super) const KIND_SNAPSHOT_FREEZE_RECEIPT: &str = "snapshot-freeze-receipt";
const TERMINAL_SCHEMA_VERSION: u32 = 1;

/// Authority domain separating this terminal kind's authority hash from raw
/// receipt identities.
const FREEZE_AUTHORITY_DOMAIN: &str = "org.frankensim.fs-session.freeze-authority.v1";

/// Armed per-session requirement state.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum SnapshotFreezeGateState {
    /// Armed after acknowledgement; activation stays refused until the
    /// matching predecessor-bound receipt is recorded.
    AwaitReceipt {
        /// Ordinal of the completed pause request this gate binds.
        request_ordinal: i64,
        /// Content hash of the acknowledged checkpoint receipt.
        acknowledgement_hash: ContentHash,
    },
    /// Receipt recorded (or exactly replayed) for the armed pause.
    Satisfied {
        /// Typed recorded receipt.
        receipt: RecordedSnapshotFreezeReceipt,
    },
}

/// The durable, validated projection of an fs-exec committed freeze receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RecordedSnapshotFreezeReceipt {
    /// Domain-separated authority over this terminal row.
    pub authority: ContentHash,
    /// Exact identity of the fs-exec committed receipt.
    pub freeze_identity: [u8; 32],
    /// Content id of the sealed solver-state envelope.
    pub sealed_content_id: [u8; 32],
    /// Transaction commitment covering identities and the drain report hash.
    pub payload_commitment: [u8; 32],
    /// Logical run whose workers drained before capture.
    pub drained_run: u64,
    /// Solver-instance generation bound at freeze time.
    pub solver_instance_generation: u64,
    /// Ordinal of the completed pause request.
    pub request_ordinal: i64,
    /// Hash of the acknowledged checkpoint receipt.
    pub acknowledgement_hash: ContentHash,
}

/// Outcome of one publication call.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SnapshotFreezePublicationDisposition {
    /// The terminal row was newly appended.
    Committed,
    /// The identical terminal row already existed; nothing changed.
    Replayed,
}

/// Returned by [`super::Governor::record_snapshot_freeze_receipt`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SnapshotFreezeReceiptWrite {
    /// Validated durable projection.
    pub recorded: RecordedSnapshotFreezeReceipt,
    /// Whether this call committed or replayed the immutable row.
    pub disposition: SnapshotFreezePublicationDisposition,
}

/// Validate that `freeze` binds to this governor, session, and completed
/// pause before any byte reaches the ledger.
pub(super) fn validate_binding(
    governor_hash: ContentHash,
    session: u64,
    freeze: &fs_exec::freeze::SnapshotFreezeReceipt,
) -> Result<(), (&'static str, SessionError)> {
    let binding = freeze.binding();
    if binding.session != session {
        return Err((
            "session",
            SessionError::SnapshotFreezeBindingMismatch {
                id: session,
                reason: "receipt names a different session",
            },
        ));
    }
    if binding.owner != *governor_hash.as_bytes() {
        return Err((
            "owner",
            SessionError::SnapshotFreezeBindingMismatch {
                id: session,
                reason: "receipt owner is not this governor identity",
            },
        ));
    }
    if freeze.disposition() != fs_exec::freeze::FreezeDisposition::BytesCommitted {
        return Err((
            "disposition",
            SessionError::SnapshotFreezeBindingMismatch {
                id: session,
                reason: "receipt does not state committed bytes",
            },
        ));
    }
    Ok(())
}

/// Canonical terminal payload bytes for one recorded receipt.
pub(super) fn encode_terminal(receipt: &RecordedSnapshotFreezeReceipt) -> Vec<u8> {
    let mut out = Vec::with_capacity(8 * 6 + 32 * 4);
    out.extend_from_slice(&TERMINAL_SCHEMA_VERSION.to_le_bytes());
    out.extend_from_slice(receipt.authority.as_bytes());
    out.extend_from_slice(&receipt.freeze_identity);
    out.extend_from_slice(&receipt.sealed_content_id);
    out.extend_from_slice(&receipt.payload_commitment);
    out.extend_from_slice(&receipt.drained_run.to_le_bytes());
    out.extend_from_slice(&receipt.solver_instance_generation.to_le_bytes());
    out.extend_from_slice(&receipt.request_ordinal.to_le_bytes());
    out.extend_from_slice(receipt.acknowledgement_hash.as_bytes());
    out
}

/// Inputs for one publication attempt; assembled under the governor lock so
/// the claim cannot drift from the armed gate.
pub(super) struct FreezePublishInputs<'a> {
    pub governor_hash: ContentHash,
    pub session_open_hash: ContentHash,
    pub session: u64,
    pub ledger_scope: &'a str,
    pub generation: u64,
    pub logical_time: i64,
    pub request_ordinal: i64,
    pub acknowledgement_hash: ContentHash,
    pub freeze: &'a fs_exec::freeze::SnapshotFreezeReceipt,
}

/// Publish the predecessor-bound receipt as ONE atomic terminal batch.
pub(super) fn publish_freeze_terminal(
    inputs: FreezePublishInputs<'_>,
    ledger: &fs_ledger::Ledger,
) -> Result<SnapshotFreezeReceiptWrite, SessionError> {
    use fs_ledger::EventRow;
    use fs_ledger::session_registry::{
        SessionMutationClaim, SessionTerminalBatch, SessionTerminalGroup, SessionTerminalRow,
    };

    if ledger.in_transaction() {
        return Err(SessionError::Persistence {
            what: "snapshot-freeze receipt refuses a caller-owned ledger transaction".to_string(),
        });
    }
    let sink = ledger
        .checked_instance_id()
        .map_err(|error| SessionError::Persistence {
            what: format!("snapshot-freeze ledger identity unavailable: {error}"),
        })?;

    let authority = fs_blake3::hash_domain(
        FREEZE_AUTHORITY_DOMAIN,
        &freeze_authority_preimage(inputs.freeze),
    );
    let recorded = RecordedSnapshotFreezeReceipt {
        authority,
        freeze_identity: *inputs.freeze.identity(),
        sealed_content_id: *inputs.freeze.sealed_content_id(),
        payload_commitment: *inputs.freeze.payload_commitment(),
        drained_run: inputs.freeze.drained_run(),
        solver_instance_generation: inputs.freeze.binding().solver_instance_generation,
        request_ordinal: inputs.request_ordinal,
        acknowledgement_hash: inputs.acknowledgement_hash,
    };
    let terminal_bytes = encode_terminal(&recorded);

    // Idempotent replay: an identical stored row means this exact publication
    // already committed; anything else under the same authority refuses.
    let existing =
        ledger
            .session_terminal(&authority)
            .map_err(|error| SessionError::Persistence {
                what: format!("snapshot-freeze terminal lookup failed: {error}"),
            })?;
    if let Some(existing) = existing {
        if existing.receipt != terminal_bytes {
            return Err(SessionError::IndeterminateMutation {
                kind: KIND_SNAPSHOT_FREEZE_RECEIPT,
                authority,
            });
        }
        return Ok(SnapshotFreezeReceiptWrite {
            recorded,
            disposition: SnapshotFreezePublicationDisposition::Replayed,
        });
    }

    let session_be = inputs.session.to_be_bytes();
    // The event payload is the ledger's JSON observability channel; the
    // authoritative receipt bytes live in the terminal row keyed by the
    // authority hash, so the event carries that key rather than a second
    // copy of the binary receipt (EventRow::payload is a JSON &str).
    let mut event_payload = String::with_capacity(96);
    event_payload.push_str("{\"schema\":\"fs-session.snapshot-freeze-event/v1\",\"authority\":\"");
    {
        use core::fmt::Write as _;
        for byte in authority.as_bytes() {
            let _ = write!(event_payload, "{byte:02x}");
        }
    }
    event_payload.push_str("\"}");
    let event = EventRow {
        session: Some(&session_be),
        t: inputs.logical_time,
        kind: KIND_SNAPSHOT_FREEZE_RECEIPT,
        payload: Some(&event_payload),
    };
    let claim = SessionMutationClaim {
        authority,
        ledger_instance_id: sink,
        governor_hash: inputs.governor_hash,
        session_open_hash: inputs.session_open_hash,
        kind: KIND_SNAPSHOT_FREEZE_RECEIPT,
        session: inputs.session,
        ledger_scope: inputs.ledger_scope,
        generation: inputs.generation,
        causal_ordinal: None,
        payload: &terminal_bytes,
    };
    let group = SessionTerminalGroup {
        terminal: SessionTerminalRow {
            claim,
            permit: None,
            receipt: &terminal_bytes,
        },
        events: &[event],
    };
    let result = ledger
        .append_session_terminal_batch(&SessionTerminalBatch { groups: &[group] })
        .map_err(|error| SessionError::Persistence {
            what: format!("snapshot-freeze terminal batch failed: {error}"),
        })?;
    match result {
        fs_ledger::session_registry::SessionTerminalBatchResult::Committed {
            terminals_inserted: 1,
            events_appended: 1,
            ..
        } => Ok(SnapshotFreezeReceiptWrite {
            recorded,
            disposition: SnapshotFreezePublicationDisposition::Committed,
        }),
        other => Err(SessionError::Persistence {
            what: format!(
                "snapshot-freeze terminal batch returned unexpected cardinality {other:?}"
            ),
        }),
    }
}

fn freeze_authority_preimage(freeze: &fs_exec::freeze::SnapshotFreezeReceipt) -> Vec<u8> {
    let mut preimage = Vec::with_capacity(128);
    preimage.extend_from_slice(freeze.identity());
    preimage.extend_from_slice(freeze.sealed_content_id());
    preimage.extend_from_slice(freeze.payload_commitment());
    preimage.extend_from_slice(&freeze.drained_run().to_le_bytes());
    preimage.extend_from_slice(&freeze.binding().solver_instance_generation.to_le_bytes());
    preimage
}

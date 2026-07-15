//! The resource GOVERNOR: continuous metering against capability tokens
//! (throttle at the grant, pause past the hard bound — NEVER a silent
//! kill), idempotency-keyed exactly-once submission, and the DECLARED
//! degradation ladder under memory pressure (spill coldest arenas →
//! coarsen adaptively → pause-serialize-resume), every event recorded
//! with attribution and flushable to the Design Ledger.

use crate::token::{CapabilityToken, SessionId};
use crate::{Guidance, SessionError};
use fs_exec::CancelGate;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

pub(crate) mod recovery;

/// Hard-bound ratio: past 6/5 of a grant the session pauses. Float and exact
/// integer resource paths derive from this one policy definition.
const HARD_FACTOR_NUMERATOR: u32 = 6;
const HARD_FACTOR_DENOMINATOR: u32 = 5;
#[allow(clippy::cast_lossless)] // small policy integers are exactly representable as f64
const HARD_FACTOR: f64 = HARD_FACTOR_NUMERATOR as f64 / HARD_FACTOR_DENOMINATOR as f64;
/// Semantic version of the restart-stable governor namespace.
pub const DURABLE_GOVERNOR_IDENTITY_VERSION: u32 = 1;
/// Semantic version of a retryable session-open authority.
pub const SESSION_OPEN_IDENTITY_VERSION: u32 = 2;
/// Semantic version of the canonical capability-token digest.
pub const SESSION_TOKEN_IDENTITY_VERSION: u32 = 1;
/// Semantic version shared by initial and resumed gate bindings.
pub const GATE_BINDING_IDENTITY_VERSION: u32 = 1;
/// Semantic version shared by client and submission meter authorities.
pub const METER_REPORT_IDENTITY_VERSION: u32 = 2;
/// Semantic version of a retryable pressure-action authority.
pub const PRESSURE_ACTION_IDENTITY_VERSION: u32 = 2;
/// Semantic version of the bounded submission-agent key digest.
pub const SUBMISSION_AGENT_KEY_IDENTITY_VERSION: u32 = 1;
/// Semantic version of the bounded canonical-program digest.
pub const SUBMISSION_PROGRAM_IDENTITY_VERSION: u32 = 1;
/// Semantic version of a retryable submission slot authority.
pub const SUBMISSION_REQUEST_IDENTITY_VERSION: u32 = 2;
/// Semantic version of retained diagnostic evidence.
pub const RETAINED_EVIDENCE_IDENTITY_VERSION: u32 = 1;
/// Semantic version of a pause-acknowledgement authority.
pub const PAUSE_ACKNOWLEDGEMENT_IDENTITY_VERSION: u32 = 1;
/// Semantic version of a resume-activation authority.
pub const RESUME_ACTIVATION_IDENTITY_VERSION: u32 = 1;
/// Semantic version of an admitted session-open receipt.
pub const SESSION_OPEN_RECEIPT_IDENTITY_VERSION: u32 = 1;
/// Semantic version of a committed meter receipt.
pub const METER_RECEIPT_IDENTITY_VERSION: u32 = 1;
/// Semantic version of a committed pressure receipt.
pub const PRESSURE_RECEIPT_IDENTITY_VERSION: u32 = 1;
/// Semantic version of a terminal submission receipt.
pub const SUBMISSION_RECEIPT_IDENTITY_VERSION: u32 = 3;
/// Semantic version of a pause-acknowledgement receipt.
pub const PAUSE_ACKNOWLEDGEMENT_RECEIPT_IDENTITY_VERSION: u32 = 1;
/// Semantic version of a resume-activation receipt.
pub const RESUME_ACTIVATION_RECEIPT_IDENTITY_VERSION: u32 = 1;
const IDEMPOTENCY_KEY_DOMAIN: &str = "org.frankensim.fs-session.idempotency-key.v3";
const IDEMPOTENCY_AGENT_DOMAIN: &str = "org.frankensim.fs-session.idempotency-agent.v1";
const IDEMPOTENCY_PROGRAM_DOMAIN: &str = "org.frankensim.fs-session.idempotency-program.v1";
const SUBMISSION_RECEIPT_DOMAIN: &str = "org.frankensim.fs-session.submission-receipt.v3";
const RETAINED_EVIDENCE_DOMAIN: &str = "org.frankensim.fs-session.retained-evidence.v1";
const SESSION_OPEN_ID_DOMAIN: &str = "org.frankensim.fs-session.open-id.v2";
const SESSION_OPEN_RECEIPT_DOMAIN: &str = "org.frankensim.fs-session.open-receipt.v1";
const SESSION_TOKEN_IDENTITY_DOMAIN: &str = "org.frankensim.fs-session.token-identity.v1";
const GATE_BINDING_ID_DOMAIN: &str = "org.frankensim.fs-session.gate-binding-id.v1";
const METER_REPORT_ID_DOMAIN: &str = "org.frankensim.fs-session.meter-report-id.v2";
const METER_RECEIPT_DOMAIN: &str = "org.frankensim.fs-session.meter-receipt.v1";
const PRESSURE_ACTION_ID_DOMAIN: &str = "org.frankensim.fs-session.pressure-action-id.v2";
const PRESSURE_RECEIPT_DOMAIN: &str = "org.frankensim.fs-session.pressure-receipt.v1";
const SUBMISSION_REQUEST_ID_DOMAIN: &str = "org.frankensim.fs-session.submission-request-id.v2";
const PAUSE_ACKNOWLEDGEMENT_ID_DOMAIN: &str =
    "org.frankensim.fs-session.pause-acknowledgement-id.v1";
const PAUSE_ACKNOWLEDGEMENT_RECEIPT_DOMAIN: &str =
    "org.frankensim.fs-session.pause-acknowledgement-receipt.v1";
const RESUME_ACTIVATION_ID_DOMAIN: &str = "org.frankensim.fs-session.resume-activation-id.v1";
const RESUME_ACTIVATION_RECEIPT_DOMAIN: &str =
    "org.frankensim.fs-session.resume-activation-receipt.v1";
const EPHEMERAL_GOVERNOR_ID_DOMAIN: &str = "org.frankensim.fs-session.ephemeral-governor-id.v1";
const DURABLE_GOVERNOR_ID_DOMAIN: &str = "org.frankensim.fs-session.durable-governor-id.v1";
const SESSION_OPEN_AUDIT_SCHEMA: &str = "fs-session-open-v1";
const SESSION_METER_AUDIT_SCHEMA: &str = "fs-session-meter-report-v1";
const SESSION_SUBMISSION_AUDIT_SCHEMA: &str = "fs-session-idempotency-v5";
const SESSION_DEGRADATION_AUDIT_SCHEMA: &str = "fs-session-degradation-v5";

// These private witness shapes make every byte-level identity input explicit
// in this owner file. The production encoders remain authoritative and are
// fingerprinted by the declarations below; the witnesses let the identity
// gate reject an unclassified field before registry data can drift.
#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct DurableGovernorIdentitySource {
    ledger_instance_id: [u8; 16],
    nonce: [u8; 32],
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct SessionOpenIdIdentitySource {
    governor_id: fs_blake3::ContentHash,
    session: u64,
    client_key_digest: fs_blake3::ContentHash,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct SessionTokenIdentitySource {
    session: u64,
    core_s_bits: u64,
    mem_bytes: u64,
    wall_s_bits: u64,
    cores: u64,
    ledger_scope: Vec<u8>,
    operations: Vec<Vec<u8>>,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
enum GateBindingIdentitySource {
    Session {
        open_id: fs_blake3::ContentHash,
    },
    Resumed {
        governor_id: fs_blake3::ContentHash,
        session: u64,
        gate_generation: u64,
        requested_ordinal: i64,
        resume_generation: u64,
    },
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
enum MeterReportIdIdentitySource {
    Client {
        governor_id: fs_blake3::ContentHash,
        session: u64,
        session_open: fs_blake3::ContentHash,
        generation: u64,
        client_key_digest: fs_blake3::ContentHash,
    },
    Submission {
        request_id: fs_blake3::ContentHash,
    },
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct PressureActionIdIdentitySource {
    governor_id: fs_blake3::ContentHash,
    session: u64,
    session_open: fs_blake3::ContentHash,
    generation: u64,
    client_key_digest: fs_blake3::ContentHash,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct SubmissionAgentKeyIdentitySource {
    value: Vec<u8>,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct SubmissionProgramIdentitySource {
    value: Vec<u8>,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct SubmissionRequestIdIdentitySource {
    governor_id: fs_blake3::ContentHash,
    session: u64,
    session_open: fs_blake3::ContentHash,
    generation: u64,
    key_hash: fs_blake3::ContentHash,
    request_hash: fs_blake3::ContentHash,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct RetainedEvidenceIdentitySource {
    value: Vec<u8>,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct PauseAcknowledgementIdIdentitySource {
    governor_id: fs_blake3::ContentHash,
    session: u64,
    gate_generation: u64,
    requested_ordinal: i64,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct ResumeActivationIdIdentitySource {
    governor_id: fs_blake3::ContentHash,
    session: u64,
    session_open: fs_blake3::ContentHash,
    acknowledgement_hash: fs_blake3::ContentHash,
    resume_generation: u64,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct SessionOpenReceiptIdentitySource {
    open_id: fs_blake3::ContentHash,
    token_digest: fs_blake3::ContentHash,
    gate_identity: Option<fs_blake3::ContentHash>,
    ledger_scope: Vec<u8>,
}

#[allow(dead_code)]
#[derive(Clone, Copy, PartialEq, Eq)]
struct ChargeIdentitySource {
    core_s_bits: u64,
    mem_peak_bytes: u64,
    wall_s_bits: u64,
}

#[allow(dead_code)]
#[derive(Clone, Copy, PartialEq, Eq)]
struct MeterSnapshotIdentitySource {
    core_s_bits: u64,
    mem_peak_bytes: u64,
    wall_s_bits: u64,
    throttled: u32,
    paused: u32,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
enum EnforcementIdentitySource {
    Ok,
    Throttled {
        resource: Vec<u8>,
        used_bits: u64,
        granted_bits: u64,
    },
    Paused {
        resource: Vec<u8>,
        used_bits: u64,
        granted_bits: u64,
        resume_hint: Vec<u8>,
    },
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct MeterReceiptIdentitySource {
    report_id: fs_blake3::ContentHash,
    commit_ordinal: u64,
    delta: ChargeIdentitySource,
    before: MeterSnapshotIdentitySource,
    after: MeterSnapshotIdentitySource,
    enforcement: EnforcementIdentitySource,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct RetainedEvidenceReceiptIdentitySource {
    preview: Vec<u8>,
    byte_len: u64,
    digest: fs_blake3::ContentHash,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct PressureEventIdentitySource {
    session: u64,
    step_tag: u8,
    pressure_level: u8,
    phase_tag: u8,
    attribution: Vec<u8>,
    ordinal: i64,
    requested_ordinal: Option<i64>,
    checkpoint: Option<RetainedEvidenceReceiptIdentitySource>,
    gate_generation: Option<u64>,
    pause_request_id: Option<PauseRequestId>,
    pressure_action_id: Option<PressureActionId>,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct PressureReceiptIdentitySource {
    action_id: fs_blake3::ContentHash,
    level: u8,
    events: Vec<PressureEventIdentitySource>,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
enum SubmissionReceiptIdentitySource {
    Done {
        request_id: fs_blake3::ContentHash,
        ledger_scope: Vec<u8>,
        admission_ordinal: u64,
        charge: ChargeIdentitySource,
        meter_receipt: fs_blake3::ContentHash,
    },
    Failed {
        request_id: fs_blake3::ContentHash,
        ledger_scope: Vec<u8>,
        admission_ordinal: u64,
        evidence: RetainedEvidenceReceiptIdentitySource,
    },
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct PauseAcknowledgementReceiptIdentitySource {
    governor_id: fs_blake3::ContentHash,
    session: u64,
    gate_generation: u64,
    requested_ordinal: i64,
    event: PressureEventIdentitySource,
    resume_generation: u64,
    gate_binding: fs_blake3::ContentHash,
}

#[allow(dead_code)]
#[derive(Clone, PartialEq, Eq)]
struct ResumeActivationReceiptIdentitySource {
    activation_id: fs_blake3::ContentHash,
    acknowledgement_hash: fs_blake3::ContentHash,
    gate_binding: fs_blake3::ContentHash,
}

#[allow(dead_code)]
fn identity_push_framed(payload: &mut Vec<u8>, bytes: &[u8]) {
    payload.extend_from_slice(
        &u64::try_from(bytes.len())
            .expect("bounded identity field fits u64")
            .to_le_bytes(),
    );
    payload.extend_from_slice(bytes);
}

#[allow(dead_code)]
fn identity_push_hash(payload: &mut Vec<u8>, hash: fs_blake3::ContentHash) {
    payload.extend_from_slice(hash.as_bytes());
}

#[allow(dead_code)]
fn identity_push_snapshot(payload: &mut Vec<u8>, snapshot: MeterSnapshotIdentitySource) {
    payload.extend_from_slice(&snapshot.core_s_bits.to_le_bytes());
    payload.extend_from_slice(&snapshot.mem_peak_bytes.to_le_bytes());
    payload.extend_from_slice(&snapshot.wall_s_bits.to_le_bytes());
    payload.extend_from_slice(&snapshot.throttled.to_le_bytes());
    payload.extend_from_slice(&snapshot.paused.to_le_bytes());
}

#[allow(dead_code)]
fn identity_push_enforcement(payload: &mut Vec<u8>, source: &EnforcementIdentitySource) {
    match source {
        EnforcementIdentitySource::Ok => payload.push(0),
        EnforcementIdentitySource::Throttled {
            resource,
            used_bits,
            granted_bits,
        } => {
            payload.push(1);
            identity_push_framed(payload, resource);
            payload.extend_from_slice(&used_bits.to_le_bytes());
            payload.extend_from_slice(&granted_bits.to_le_bytes());
        }
        EnforcementIdentitySource::Paused {
            resource,
            used_bits,
            granted_bits,
            resume_hint,
        } => {
            payload.push(2);
            identity_push_framed(payload, resource);
            payload.extend_from_slice(&used_bits.to_le_bytes());
            payload.extend_from_slice(&granted_bits.to_le_bytes());
            identity_push_framed(payload, resume_hint);
        }
    }
}

#[allow(dead_code)]
fn identity_push_pressure_event(payload: &mut Vec<u8>, source: &PressureEventIdentitySource) {
    payload.extend_from_slice(&source.session.to_le_bytes());
    payload.push(source.step_tag);
    payload.push(source.pressure_level);
    payload.push(source.phase_tag);
    identity_push_framed(payload, &source.attribution);
    payload.extend_from_slice(&source.ordinal.to_le_bytes());
    match source.requested_ordinal {
        Some(value) => {
            payload.push(1);
            payload.extend_from_slice(&value.to_le_bytes());
        }
        None => payload.push(0),
    }
    match &source.checkpoint {
        Some(checkpoint) => {
            payload.push(1);
            payload.extend_from_slice(&checkpoint.byte_len.to_le_bytes());
            identity_push_hash(payload, checkpoint.digest);
        }
        None => payload.push(0),
    }
    match source.gate_generation {
        Some(value) => {
            payload.push(1);
            payload.extend_from_slice(&value.to_le_bytes());
        }
        None => payload.push(0),
    }
}

#[allow(dead_code)]
fn durable_governor_identity(source: &DurableGovernorIdentitySource) -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    payload.extend_from_slice(&source.ledger_instance_id);
    payload.extend_from_slice(&source.nonce);
    fs_blake3::hash_domain(DURABLE_GOVERNOR_ID_DOMAIN, &payload)
}

#[allow(dead_code)]
fn session_open_authority_identity(source: &SessionOpenIdIdentitySource) -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    identity_push_hash(&mut payload, source.governor_id);
    payload.extend_from_slice(&source.session.to_le_bytes());
    identity_push_hash(&mut payload, source.client_key_digest);
    fs_blake3::hash_domain(SESSION_OPEN_ID_DOMAIN, &payload)
}

#[allow(dead_code)]
fn session_token_identity(source: &SessionTokenIdentitySource) -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    payload.extend_from_slice(&source.session.to_le_bytes());
    payload.extend_from_slice(&source.core_s_bits.to_le_bytes());
    payload.extend_from_slice(&source.mem_bytes.to_le_bytes());
    payload.extend_from_slice(&source.wall_s_bits.to_le_bytes());
    payload.extend_from_slice(&source.cores.to_le_bytes());
    identity_push_framed(&mut payload, &source.ledger_scope);
    payload.extend_from_slice(
        &u64::try_from(source.operations.len())
            .expect("bounded operator count fits u64")
            .to_le_bytes(),
    );
    for operation in &source.operations {
        identity_push_framed(&mut payload, operation);
    }
    fs_blake3::hash_domain(SESSION_TOKEN_IDENTITY_DOMAIN, &payload)
}

#[allow(dead_code)]
fn gate_binding_identity(source: &GateBindingIdentitySource) -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    match source {
        GateBindingIdentitySource::Session { open_id } => {
            identity_push_hash(&mut payload, *open_id);
        }
        GateBindingIdentitySource::Resumed {
            governor_id,
            session,
            gate_generation,
            requested_ordinal,
            resume_generation,
        } => {
            identity_push_hash(&mut payload, *governor_id);
            payload.extend_from_slice(&session.to_le_bytes());
            payload.extend_from_slice(&gate_generation.to_le_bytes());
            payload.extend_from_slice(&requested_ordinal.to_le_bytes());
            payload.extend_from_slice(&resume_generation.to_le_bytes());
        }
    }
    fs_blake3::hash_domain(GATE_BINDING_ID_DOMAIN, &payload)
}

#[allow(dead_code)]
fn meter_report_authority_identity(source: &MeterReportIdIdentitySource) -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    match source {
        MeterReportIdIdentitySource::Client {
            governor_id,
            session,
            session_open,
            generation,
            client_key_digest,
        } => {
            identity_push_hash(&mut payload, *governor_id);
            payload.extend_from_slice(&session.to_le_bytes());
            identity_push_hash(&mut payload, *session_open);
            payload.extend_from_slice(&generation.to_le_bytes());
            identity_push_hash(&mut payload, *client_key_digest);
        }
        MeterReportIdIdentitySource::Submission { request_id } => {
            identity_push_hash(&mut payload, *request_id);
        }
    }
    fs_blake3::hash_domain(METER_REPORT_ID_DOMAIN, &payload)
}

#[allow(dead_code)]
fn pressure_action_authority_identity(
    source: &PressureActionIdIdentitySource,
) -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    identity_push_hash(&mut payload, source.governor_id);
    payload.extend_from_slice(&source.session.to_le_bytes());
    identity_push_hash(&mut payload, source.session_open);
    payload.extend_from_slice(&source.generation.to_le_bytes());
    identity_push_hash(&mut payload, source.client_key_digest);
    fs_blake3::hash_domain(PRESSURE_ACTION_ID_DOMAIN, &payload)
}

#[allow(dead_code)]
fn submission_agent_key_identity(
    source: &SubmissionAgentKeyIdentitySource,
) -> fs_blake3::ContentHash {
    fs_blake3::hash_domain(IDEMPOTENCY_AGENT_DOMAIN, &source.value)
}

#[allow(dead_code)]
fn submission_program_identity(source: &SubmissionProgramIdentitySource) -> fs_blake3::ContentHash {
    fs_blake3::hash_domain(IDEMPOTENCY_PROGRAM_DOMAIN, &source.value)
}

#[allow(dead_code)]
fn submission_request_authority_identity(
    source: &SubmissionRequestIdIdentitySource,
) -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    identity_push_hash(&mut payload, source.governor_id);
    payload.extend_from_slice(&source.session.to_le_bytes());
    identity_push_hash(&mut payload, source.session_open);
    payload.extend_from_slice(&source.generation.to_le_bytes());
    identity_push_hash(&mut payload, source.key_hash);
    fs_blake3::hash_domain(SUBMISSION_REQUEST_ID_DOMAIN, &payload)
}

#[allow(dead_code)]
fn retained_evidence_identity(source: &RetainedEvidenceIdentitySource) -> fs_blake3::ContentHash {
    fs_blake3::hash_domain(RETAINED_EVIDENCE_DOMAIN, &source.value)
}

#[allow(dead_code)]
fn pause_acknowledgement_authority_identity(
    source: &PauseAcknowledgementIdIdentitySource,
) -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    identity_push_hash(&mut payload, source.governor_id);
    payload.extend_from_slice(&source.session.to_le_bytes());
    payload.extend_from_slice(&source.gate_generation.to_le_bytes());
    payload.extend_from_slice(&source.requested_ordinal.to_le_bytes());
    fs_blake3::hash_domain(PAUSE_ACKNOWLEDGEMENT_ID_DOMAIN, &payload)
}

#[allow(dead_code)]
fn resume_activation_authority_identity(
    source: &ResumeActivationIdIdentitySource,
) -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    identity_push_hash(&mut payload, source.governor_id);
    payload.extend_from_slice(&source.session.to_le_bytes());
    identity_push_hash(&mut payload, source.session_open);
    identity_push_hash(&mut payload, source.acknowledgement_hash);
    payload.extend_from_slice(&source.resume_generation.to_le_bytes());
    fs_blake3::hash_domain(RESUME_ACTIVATION_ID_DOMAIN, &payload)
}

#[allow(dead_code)]
fn session_open_receipt_identity(
    source: &SessionOpenReceiptIdentitySource,
) -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    identity_push_hash(&mut payload, source.open_id);
    identity_push_hash(&mut payload, source.token_digest);
    match source.gate_identity {
        Some(identity) => {
            payload.push(1);
            identity_push_hash(&mut payload, identity);
        }
        None => payload.push(0),
    }
    identity_push_framed(&mut payload, &source.ledger_scope);
    fs_blake3::hash_domain(SESSION_OPEN_RECEIPT_DOMAIN, &payload)
}

#[allow(dead_code)]
fn meter_receipt_identity(source: &MeterReceiptIdentitySource) -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    identity_push_hash(&mut payload, source.report_id);
    payload.extend_from_slice(&source.commit_ordinal.to_le_bytes());
    payload.extend_from_slice(&source.delta.core_s_bits.to_le_bytes());
    payload.extend_from_slice(&source.delta.mem_peak_bytes.to_le_bytes());
    payload.extend_from_slice(&source.delta.wall_s_bits.to_le_bytes());
    identity_push_snapshot(&mut payload, source.before);
    identity_push_snapshot(&mut payload, source.after);
    identity_push_enforcement(&mut payload, &source.enforcement);
    fs_blake3::hash_domain(METER_RECEIPT_DOMAIN, &payload)
}

#[allow(dead_code)]
fn pressure_receipt_identity(source: &PressureReceiptIdentitySource) -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    identity_push_hash(&mut payload, source.action_id);
    payload.push(source.level);
    payload.extend_from_slice(
        &u64::try_from(source.events.len())
            .expect("bounded degradation event count fits u64")
            .to_le_bytes(),
    );
    for event in &source.events {
        identity_push_pressure_event(&mut payload, event);
    }
    fs_blake3::hash_domain(PRESSURE_RECEIPT_DOMAIN, &payload)
}

#[allow(dead_code)]
fn submission_receipt_identity(source: &SubmissionReceiptIdentitySource) -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    match source {
        SubmissionReceiptIdentitySource::Done {
            request_id,
            ledger_scope,
            admission_ordinal,
            charge,
            meter_receipt,
        } => {
            identity_push_hash(&mut payload, *request_id);
            identity_push_framed(&mut payload, ledger_scope);
            payload.extend_from_slice(&admission_ordinal.to_le_bytes());
            payload.push(0);
            payload.extend_from_slice(&charge.core_s_bits.to_le_bytes());
            payload.extend_from_slice(&charge.mem_peak_bytes.to_le_bytes());
            payload.extend_from_slice(&charge.wall_s_bits.to_le_bytes());
            identity_push_hash(&mut payload, *meter_receipt);
        }
        SubmissionReceiptIdentitySource::Failed {
            request_id,
            ledger_scope,
            admission_ordinal,
            evidence,
        } => {
            identity_push_hash(&mut payload, *request_id);
            identity_push_framed(&mut payload, ledger_scope);
            payload.extend_from_slice(&admission_ordinal.to_le_bytes());
            payload.push(1);
            payload.extend_from_slice(&evidence.byte_len.to_le_bytes());
            identity_push_hash(&mut payload, evidence.digest);
        }
    }
    fs_blake3::hash_domain(SUBMISSION_RECEIPT_DOMAIN, &payload)
}

#[allow(dead_code)]
fn pause_acknowledgement_receipt_identity(
    source: &PauseAcknowledgementReceiptIdentitySource,
) -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    identity_push_hash(&mut payload, source.governor_id);
    payload.extend_from_slice(&source.session.to_le_bytes());
    payload.extend_from_slice(&source.gate_generation.to_le_bytes());
    payload.extend_from_slice(&source.requested_ordinal.to_le_bytes());
    identity_push_pressure_event(&mut payload, &source.event);
    payload.extend_from_slice(&source.resume_generation.to_le_bytes());
    identity_push_hash(&mut payload, source.gate_binding);
    fs_blake3::hash_domain(PAUSE_ACKNOWLEDGEMENT_RECEIPT_DOMAIN, &payload)
}

#[allow(dead_code)]
fn resume_activation_receipt_identity(
    source: &ResumeActivationReceiptIdentitySource,
) -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    identity_push_hash(&mut payload, source.activation_id);
    identity_push_hash(&mut payload, source.acknowledgement_hash);
    identity_push_hash(&mut payload, source.gate_binding);
    fs_blake3::hash_domain(RESUME_ACTIVATION_RECEIPT_DOMAIN, &payload)
}

#[allow(dead_code)]
fn session_identity_transport_is_current(
    terminal_schema_version: u32,
    audit_schema: Option<&str>,
) -> bool {
    if terminal_schema_version != recovery::TERMINAL_SCHEMA_VERSION {
        return false;
    }
    match audit_schema {
        None => true,
        Some(schema) => matches!(
            schema,
            SESSION_OPEN_AUDIT_SCHEMA
                | SESSION_METER_AUDIT_SCHEMA
                | SESSION_SUBMISSION_AUDIT_SCHEMA
                | SESSION_DEGRADATION_AUDIT_SCHEMA
        ),
    }
}

#[allow(dead_code)]
fn classify_durable_governor_identity_fields(source: &DurableGovernorIdentitySource) {
    let DurableGovernorIdentitySource {
        ledger_instance_id,
        nonce,
    } = source;
    let _ = (ledger_instance_id, nonce);
}

#[allow(dead_code)]
fn classify_session_open_id_identity_fields(source: &SessionOpenIdIdentitySource) {
    let SessionOpenIdIdentitySource {
        governor_id,
        session,
        client_key_digest,
    } = source;
    let _ = (governor_id, session, client_key_digest);
}

#[allow(dead_code)]
fn classify_session_token_identity_fields(source: &SessionTokenIdentitySource) {
    let SessionTokenIdentitySource {
        session,
        core_s_bits,
        mem_bytes,
        wall_s_bits,
        cores,
        ledger_scope,
        operations,
    } = source;
    let _ = (
        session,
        core_s_bits,
        mem_bytes,
        wall_s_bits,
        cores,
        ledger_scope,
        operations,
    );
}

#[allow(dead_code)]
fn classify_gate_binding_identity_fields(source: &GateBindingIdentitySource) {
    match source {
        GateBindingIdentitySource::Session { open_id } => {
            let _ = open_id;
        }
        GateBindingIdentitySource::Resumed {
            governor_id,
            session,
            gate_generation,
            requested_ordinal,
            resume_generation,
        } => {
            let _ = (
                governor_id,
                session,
                gate_generation,
                requested_ordinal,
                resume_generation,
            );
        }
    }
}

#[allow(dead_code)]
fn classify_meter_report_id_identity_fields(source: &MeterReportIdIdentitySource) {
    match source {
        MeterReportIdIdentitySource::Client {
            governor_id,
            session,
            session_open,
            generation,
            client_key_digest,
        } => {
            let _ = (
                governor_id,
                session,
                session_open,
                generation,
                client_key_digest,
            );
        }
        MeterReportIdIdentitySource::Submission { request_id } => {
            let _ = request_id;
        }
    }
}

#[allow(dead_code)]
fn classify_pressure_action_id_identity_fields(source: &PressureActionIdIdentitySource) {
    let PressureActionIdIdentitySource {
        governor_id,
        session,
        session_open,
        generation,
        client_key_digest,
    } = source;
    let _ = (
        governor_id,
        session,
        session_open,
        generation,
        client_key_digest,
    );
}

#[allow(dead_code)]
fn classify_submission_agent_key_identity_fields(source: &SubmissionAgentKeyIdentitySource) {
    let SubmissionAgentKeyIdentitySource { value } = source;
    let _ = value;
}

#[allow(dead_code)]
fn classify_submission_program_identity_fields(source: &SubmissionProgramIdentitySource) {
    let SubmissionProgramIdentitySource { value } = source;
    let _ = value;
}

#[allow(dead_code)]
fn classify_submission_request_id_identity_fields(source: &SubmissionRequestIdIdentitySource) {
    let SubmissionRequestIdIdentitySource {
        governor_id,
        session,
        session_open,
        generation,
        key_hash,
        request_hash,
    } = source;
    let _ = (
        governor_id,
        session,
        session_open,
        generation,
        key_hash,
        request_hash,
    );
}

#[allow(dead_code)]
fn classify_retained_evidence_identity_fields(source: &RetainedEvidenceIdentitySource) {
    let RetainedEvidenceIdentitySource { value } = source;
    let _ = value;
}

#[allow(dead_code)]
fn classify_pause_acknowledgement_id_identity_fields(
    source: &PauseAcknowledgementIdIdentitySource,
) {
    let PauseAcknowledgementIdIdentitySource {
        governor_id,
        session,
        gate_generation,
        requested_ordinal,
    } = source;
    let _ = (governor_id, session, gate_generation, requested_ordinal);
}

#[allow(dead_code)]
fn classify_resume_activation_id_identity_fields(source: &ResumeActivationIdIdentitySource) {
    let ResumeActivationIdIdentitySource {
        governor_id,
        session,
        session_open,
        acknowledgement_hash,
        resume_generation,
    } = source;
    let _ = (
        governor_id,
        session,
        session_open,
        acknowledgement_hash,
        resume_generation,
    );
}

#[allow(dead_code)]
fn classify_session_open_receipt_identity_fields(source: &SessionOpenReceiptIdentitySource) {
    let SessionOpenReceiptIdentitySource {
        open_id,
        token_digest,
        gate_identity,
        ledger_scope,
    } = source;
    let _ = (open_id, token_digest, gate_identity, ledger_scope);
}

#[allow(dead_code)]
fn classify_meter_receipt_identity_fields(
    source: &MeterReceiptIdentitySource,
    charge: &ChargeIdentitySource,
    snapshot: &MeterSnapshotIdentitySource,
    enforcement: &EnforcementIdentitySource,
) {
    let MeterReceiptIdentitySource {
        report_id,
        commit_ordinal,
        delta,
        before,
        after,
        enforcement: receipt_enforcement,
    } = source;
    let ChargeIdentitySource {
        core_s_bits,
        mem_peak_bytes,
        wall_s_bits,
    } = charge;
    let MeterSnapshotIdentitySource {
        core_s_bits: snapshot_core_s_bits,
        mem_peak_bytes: snapshot_mem_peak_bytes,
        wall_s_bits: snapshot_wall_s_bits,
        throttled,
        paused,
    } = snapshot;
    match enforcement {
        EnforcementIdentitySource::Ok => {}
        EnforcementIdentitySource::Throttled {
            resource,
            used_bits,
            granted_bits,
        } => {
            let _ = (resource, used_bits, granted_bits);
        }
        EnforcementIdentitySource::Paused {
            resource,
            used_bits,
            granted_bits,
            resume_hint,
        } => {
            let _ = (resource, used_bits, granted_bits, resume_hint);
        }
    }
    let _ = (
        report_id,
        commit_ordinal,
        delta,
        before,
        after,
        receipt_enforcement,
        core_s_bits,
        mem_peak_bytes,
        wall_s_bits,
        snapshot_core_s_bits,
        snapshot_mem_peak_bytes,
        snapshot_wall_s_bits,
        throttled,
        paused,
    );
}

#[allow(dead_code)]
fn classify_retained_evidence_receipt_identity_fields(
    source: &RetainedEvidenceReceiptIdentitySource,
) {
    let RetainedEvidenceReceiptIdentitySource {
        preview,
        byte_len,
        digest,
    } = source;
    let _ = (preview, byte_len, digest);
}

#[allow(dead_code)]
fn classify_pressure_event_identity_fields(source: &PressureEventIdentitySource) {
    let PressureEventIdentitySource {
        session,
        step_tag,
        pressure_level,
        phase_tag,
        attribution,
        ordinal,
        requested_ordinal,
        checkpoint,
        gate_generation,
        pause_request_id,
        pressure_action_id,
    } = source;
    let _ = (
        session,
        step_tag,
        pressure_level,
        phase_tag,
        attribution,
        ordinal,
        requested_ordinal,
        checkpoint,
        gate_generation,
        pause_request_id,
        pressure_action_id,
    );
}

#[allow(dead_code)]
fn classify_pressure_receipt_identity_fields(
    source: &PressureReceiptIdentitySource,
    event: &PressureEventIdentitySource,
    evidence: &RetainedEvidenceReceiptIdentitySource,
) {
    let PressureReceiptIdentitySource {
        action_id,
        level,
        events,
    } = source;
    let PressureEventIdentitySource {
        session,
        step_tag,
        pressure_level,
        phase_tag,
        attribution,
        ordinal,
        requested_ordinal,
        checkpoint,
        gate_generation,
        pause_request_id,
        pressure_action_id,
    } = event;
    let RetainedEvidenceReceiptIdentitySource {
        preview,
        byte_len,
        digest,
    } = evidence;
    let _ = (
        action_id,
        level,
        events,
        session,
        step_tag,
        pressure_level,
        phase_tag,
        attribution,
        ordinal,
        requested_ordinal,
        checkpoint,
        gate_generation,
        pause_request_id,
        pressure_action_id,
        preview,
        byte_len,
        digest,
    );
}

#[allow(dead_code)]
fn classify_submission_receipt_identity_fields(
    source: &SubmissionReceiptIdentitySource,
    charge: &ChargeIdentitySource,
    evidence: &RetainedEvidenceReceiptIdentitySource,
) {
    match source {
        SubmissionReceiptIdentitySource::Done {
            request_id,
            ledger_scope,
            admission_ordinal,
            charge: terminal_charge,
            meter_receipt,
        } => {
            let _ = (
                request_id,
                ledger_scope,
                admission_ordinal,
                terminal_charge,
                meter_receipt,
            );
        }
        SubmissionReceiptIdentitySource::Failed {
            request_id,
            ledger_scope,
            admission_ordinal,
            evidence: terminal_evidence,
        } => {
            let _ = (
                request_id,
                ledger_scope,
                admission_ordinal,
                terminal_evidence,
            );
        }
    }
    let ChargeIdentitySource {
        core_s_bits,
        mem_peak_bytes,
        wall_s_bits,
    } = charge;
    let RetainedEvidenceReceiptIdentitySource {
        preview,
        byte_len,
        digest,
    } = evidence;
    let _ = (
        core_s_bits,
        mem_peak_bytes,
        wall_s_bits,
        preview,
        byte_len,
        digest,
    );
}

#[allow(dead_code)]
fn classify_pause_acknowledgement_receipt_identity_fields(
    source: &PauseAcknowledgementReceiptIdentitySource,
    event: &PressureEventIdentitySource,
    evidence: &RetainedEvidenceReceiptIdentitySource,
) {
    let PauseAcknowledgementReceiptIdentitySource {
        governor_id,
        session,
        gate_generation,
        requested_ordinal,
        event: receipt_event,
        resume_generation,
        gate_binding,
    } = source;
    let PressureEventIdentitySource {
        session: event_session,
        step_tag,
        pressure_level,
        phase_tag,
        attribution,
        ordinal,
        requested_ordinal: event_requested_ordinal,
        checkpoint,
        gate_generation: event_gate_generation,
        pause_request_id,
        pressure_action_id,
    } = event;
    let RetainedEvidenceReceiptIdentitySource {
        preview,
        byte_len,
        digest,
    } = evidence;
    let _ = (
        governor_id,
        session,
        gate_generation,
        requested_ordinal,
        receipt_event,
        resume_generation,
        gate_binding,
        event_session,
        step_tag,
        pressure_level,
        phase_tag,
        attribution,
        ordinal,
        event_requested_ordinal,
        checkpoint,
        event_gate_generation,
        pause_request_id,
        pressure_action_id,
        preview,
        byte_len,
        digest,
    );
}

#[allow(dead_code)]
fn classify_resume_activation_receipt_identity_fields(
    source: &ResumeActivationReceiptIdentitySource,
) {
    let ResumeActivationReceiptIdentitySource {
        activation_id,
        acknowledgement_hash,
        gate_binding,
    } = source;
    let _ = (activation_id, acknowledgement_hash, gate_binding);
}

/// Owner-local durable-governor declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const DURABLE_GOVERNOR_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-session:durable-governor-id",
    "version_const=DURABLE_GOVERNOR_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-session.durable-governor-id.v1",
    "domain_const=DURABLE_GOVERNOR_ID_DOMAIN",
    "encoder=durable_governor_identity",
    "encoder_helpers=none",
    "schema_constants=DURABLE_GOVERNOR_IDENTITY_VERSION,DURABLE_GOVERNOR_ID_DOMAIN,crates/fs-session/src/governor/recovery.rs#TERMINAL_SCHEMA_VERSION",
    "schema_functions=Governor::new_durable,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=fs-ledger:physical-instance",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=DurableGovernorIdentitySource",
    "source_fields=DurableGovernorIdentitySource.ledger_instance_id:semantic,DurableGovernorIdentitySource.nonce:semantic",
    "source_bindings=DurableGovernorIdentitySource.ledger_instance_id>ledger-instance-id,DurableGovernorIdentitySource.nonce>durable-nonce",
    "external_semantic_fields=identity-version,digest-domain,canonical-field-order",
    "semantic_fields=identity-version,digest-domain,canonical-field-order,ledger-instance-id,durable-nonce",
    "excluded_fields=none",
    "consumers=Governor::identity,Governor::session_open_id,Governor::meter_report_id,Governor::submission_request_id,Governor::pressure_action_id",
    "mutations=identity-version:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,digest-domain:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,canonical-field-order:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,ledger-instance-id:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,durable-nonce:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_durable_governor_identity_fields",
    "transport_guard=session_identity_transport_is_current",
    "version_guard=crates/fs-session/src/governor.rs#governor_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-session:durable-governor-id",
];

/// Owner-local session-open authority declaration.
#[allow(dead_code)]
pub const SESSION_OPEN_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-session:session-open-id",
    "version_const=SESSION_OPEN_IDENTITY_VERSION",
    "version=2",
    "domain=org.frankensim.fs-session.open-id.v2",
    "domain_const=SESSION_OPEN_ID_DOMAIN",
    "encoder=session_open_authority_identity",
    "encoder_helpers=identity_push_hash",
    "schema_constants=SESSION_OPEN_IDENTITY_VERSION,SESSION_OPEN_ID_DOMAIN,MAX_IDEMPOTENCY_INPUT_BYTES,crates/fs-session/src/governor/recovery.rs#TERMINAL_SCHEMA_VERSION",
    "schema_functions=Governor::session_open_id,bounded_request_digest,crates/fs-session/src/governor/recovery.rs#encode_open_payload,crates/fs-session/src/governor/recovery.rs#decode_open_payload,crates/fs-session/src/governor/recovery.rs#encode_open_receipt,crates/fs-session/src/governor/recovery.rs#decode_open_receipt,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=fs-session:durable-governor-id",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=SessionOpenIdIdentitySource",
    "source_fields=SessionOpenIdIdentitySource.governor_id:semantic,SessionOpenIdIdentitySource.session:semantic,SessionOpenIdIdentitySource.client_key_digest:semantic",
    "source_bindings=SessionOpenIdIdentitySource.governor_id>governor-id,SessionOpenIdIdentitySource.session>session,SessionOpenIdIdentitySource.client_key_digest>client-key-digest",
    "external_semantic_fields=identity-version,digest-domain,canonical-field-order",
    "semantic_fields=identity-version,digest-domain,canonical-field-order,governor-id,session,client-key-digest",
    "excluded_fields=none",
    "consumers=Governor::open_session,Governor::open_session_gated,SessionOpenReceipt,fs-ledger::session_registry::SessionMutationClaim",
    "mutations=identity-version:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,digest-domain:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,canonical-field-order:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,governor-id:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,session:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,client-key-digest:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_session_open_id_identity_fields",
    "transport_guard=session_identity_transport_is_current",
    "version_guard=crates/fs-session/src/governor.rs#governor_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-session:session-open-id",
];

/// Owner-local capability-token identity declaration.
#[allow(dead_code)]
pub const SESSION_TOKEN_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-session:session-token-identity",
    "version_const=SESSION_TOKEN_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-session.token-identity.v1",
    "domain_const=SESSION_TOKEN_IDENTITY_DOMAIN",
    "encoder=session_token_identity",
    "encoder_helpers=identity_push_framed",
    "schema_constants=SESSION_TOKEN_IDENTITY_VERSION,SESSION_TOKEN_IDENTITY_DOMAIN,crates/fs-session/src/governor/recovery.rs#TERMINAL_SCHEMA_VERSION",
    "schema_functions=capability_token_identity,crates/fs-session/src/token.rs#CapabilityToken::validate_operator_grants,crates/fs-session/src/token.rs#CapabilityToken::validate_ledger_scope,crates/fs-session/src/governor/recovery.rs#encode_open_payload,crates/fs-session/src/governor/recovery.rs#decode_open_payload,buffered_open_receipt,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=none",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=SessionTokenIdentitySource",
    "source_fields=SessionTokenIdentitySource.session:semantic,SessionTokenIdentitySource.core_s_bits:semantic,SessionTokenIdentitySource.mem_bytes:semantic,SessionTokenIdentitySource.wall_s_bits:semantic,SessionTokenIdentitySource.cores:semantic,SessionTokenIdentitySource.ledger_scope:semantic,SessionTokenIdentitySource.operations:semantic",
    "source_bindings=SessionTokenIdentitySource.session>session,SessionTokenIdentitySource.core_s_bits>core-seconds-bits,SessionTokenIdentitySource.mem_bytes>memory-bytes,SessionTokenIdentitySource.wall_s_bits>wall-seconds-bits,SessionTokenIdentitySource.cores>cores,SessionTokenIdentitySource.ledger_scope>ledger-scope,SessionTokenIdentitySource.operations>operation-count+operation-order+operation-bytes",
    "external_semantic_fields=identity-version,digest-domain,canonical-field-order,length-framing",
    "semantic_fields=identity-version,digest-domain,canonical-field-order,length-framing,session,core-seconds-bits,memory-bytes,wall-seconds-bits,cores,ledger-scope,operation-count,operation-order,operation-bytes",
    "excluded_fields=none",
    "consumers=Governor::open_session,Governor::open_session_gated,SessionOpenReceipt::token_digest,buffered_open_receipt",
    "mutations=identity-version:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,digest-domain:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,canonical-field-order:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,length-framing:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,session:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,core-seconds-bits:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,memory-bytes:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,wall-seconds-bits:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,cores:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,ledger-scope:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,operation-count:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,operation-order:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,operation-bytes:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_session_token_identity_fields",
    "transport_guard=session_identity_transport_is_current",
    "version_guard=crates/fs-session/src/governor.rs#governor_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-session:session-token-identity",
];

/// Owner-local gate-binding identity declaration.
#[allow(dead_code)]
pub const GATE_BINDING_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-session:gate-binding-id",
    "version_const=GATE_BINDING_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-session.gate-binding-id.v1",
    "domain_const=GATE_BINDING_ID_DOMAIN",
    "encoder=gate_binding_identity",
    "encoder_helpers=identity_push_hash",
    "schema_constants=GATE_BINDING_IDENTITY_VERSION,GATE_BINDING_ID_DOMAIN,crates/fs-session/src/governor/recovery.rs#TERMINAL_SCHEMA_VERSION",
    "schema_functions=session_gate_binding,resumed_gate_binding,crates/fs-session/src/governor/recovery.rs#encode_open_payload,crates/fs-session/src/governor/recovery.rs#decode_open_payload,crates/fs-session/src/governor/recovery.rs#encode_pause_ack_terminal_receipt,crates/fs-session/src/governor/recovery.rs#decode_pause_ack_terminal_receipt,crates/fs-session/src/governor/recovery.rs#encode_activation_payload,crates/fs-session/src/governor/recovery.rs#decode_activation_payload,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=fs-session:pause-acknowledgement-id,fs-session:session-open-id",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=GateBindingIdentitySource",
    "source_fields=GateBindingIdentitySource.variant:semantic,GateBindingIdentitySource.open_id:semantic,GateBindingIdentitySource.governor_id:semantic,GateBindingIdentitySource.session:semantic,GateBindingIdentitySource.gate_generation:semantic,GateBindingIdentitySource.requested_ordinal:semantic,GateBindingIdentitySource.resume_generation:semantic",
    "source_bindings=GateBindingIdentitySource.variant>shape-tag,GateBindingIdentitySource.open_id>open-id,GateBindingIdentitySource.governor_id>governor-id,GateBindingIdentitySource.session>session,GateBindingIdentitySource.gate_generation>gate-generation,GateBindingIdentitySource.requested_ordinal>requested-ordinal,GateBindingIdentitySource.resume_generation>resume-generation",
    "external_semantic_fields=identity-version,digest-domain,canonical-field-order",
    "semantic_fields=identity-version,digest-domain,canonical-field-order,shape-tag,open-id,governor-id,session,gate-generation,requested-ordinal,resume-generation",
    "excluded_fields=none",
    "consumers=SessionOpenReceipt::gate_identity,PauseAcknowledgement::gate_binding,ResumeActivationReceipt::gate_binding",
    "mutations=identity-version:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,digest-domain:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,canonical-field-order:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,shape-tag:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,open-id:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,governor-id:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,session:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,gate-generation:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,requested-ordinal:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,resume-generation:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_gate_binding_identity_fields",
    "transport_guard=session_identity_transport_is_current",
    "version_guard=crates/fs-session/src/governor.rs#governor_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-session:gate-binding-id",
];

/// Owner-local meter-report authority declaration.
#[allow(dead_code)]
pub const METER_REPORT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-session:meter-report-id",
    "version_const=METER_REPORT_IDENTITY_VERSION",
    "version=2",
    "domain=org.frankensim.fs-session.meter-report-id.v2",
    "domain_const=METER_REPORT_ID_DOMAIN",
    "encoder=meter_report_authority_identity",
    "encoder_helpers=identity_push_hash",
    "schema_constants=METER_REPORT_IDENTITY_VERSION,METER_REPORT_ID_DOMAIN,MAX_IDEMPOTENCY_INPUT_BYTES,crates/fs-session/src/governor/recovery.rs#TERMINAL_SCHEMA_VERSION",
    "schema_functions=Governor::meter_report_id,Governor::submission_meter_report_id,bounded_request_digest,crates/fs-session/src/governor/recovery.rs#encode_meter_payload,crates/fs-session/src/governor/recovery.rs#decode_meter_payload,crates/fs-session/src/governor/recovery.rs#encode_meter_terminal_receipt,crates/fs-session/src/governor/recovery.rs#decode_meter_terminal_receipt,buffered_meter_receipt,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=fs-session:session-open-id,fs-session:submission-request-id",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=MeterReportIdIdentitySource",
    "source_fields=MeterReportIdIdentitySource.variant:semantic,MeterReportIdIdentitySource.governor_id:semantic,MeterReportIdIdentitySource.session:semantic,MeterReportIdIdentitySource.session_open:semantic,MeterReportIdIdentitySource.generation:semantic,MeterReportIdIdentitySource.client_key_digest:semantic,MeterReportIdIdentitySource.request_id:semantic",
    "source_bindings=MeterReportIdIdentitySource.variant>shape-tag,MeterReportIdIdentitySource.governor_id>governor-id,MeterReportIdIdentitySource.session>session,MeterReportIdIdentitySource.session_open>session-open,MeterReportIdIdentitySource.generation>generation,MeterReportIdIdentitySource.client_key_digest>client-key-digest,MeterReportIdIdentitySource.request_id>submission-request-id",
    "external_semantic_fields=identity-version,digest-domain,canonical-field-order",
    "semantic_fields=identity-version,digest-domain,canonical-field-order,shape-tag,governor-id,session,session-open,generation,client-key-digest,submission-request-id",
    "excluded_fields=none",
    "consumers=Governor::charge,MeterReceipt::report_id,buffered_meter_receipt,SubmissionReceipt",
    "mutations=identity-version:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,digest-domain:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,canonical-field-order:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,shape-tag:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,governor-id:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,session:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,session-open:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,generation:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,client-key-digest:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,submission-request-id:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_meter_report_id_identity_fields",
    "transport_guard=session_identity_transport_is_current",
    "version_guard=crates/fs-session/src/governor.rs#governor_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-session:meter-report-id",
];

/// Owner-local pressure-action authority declaration.
#[allow(dead_code)]
pub const PRESSURE_ACTION_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-session:pressure-action-id",
    "version_const=PRESSURE_ACTION_IDENTITY_VERSION",
    "version=2",
    "domain=org.frankensim.fs-session.pressure-action-id.v2",
    "domain_const=PRESSURE_ACTION_ID_DOMAIN",
    "encoder=pressure_action_authority_identity",
    "encoder_helpers=identity_push_hash",
    "schema_constants=PRESSURE_ACTION_IDENTITY_VERSION,PRESSURE_ACTION_ID_DOMAIN,MAX_IDEMPOTENCY_INPUT_BYTES,crates/fs-session/src/governor/recovery.rs#TERMINAL_SCHEMA_VERSION",
    "schema_functions=Governor::pressure_action_id,bounded_request_digest,crates/fs-session/src/governor/recovery.rs#encode_pressure_payload,crates/fs-session/src/governor/recovery.rs#decode_pressure_payload,crates/fs-session/src/governor/recovery.rs#encode_pressure_terminal_receipt,crates/fs-session/src/governor/recovery.rs#decode_pressure_terminal_receipt,buffered_degradation_event,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=fs-session:session-open-id",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=PressureActionIdIdentitySource",
    "source_fields=PressureActionIdIdentitySource.governor_id:semantic,PressureActionIdIdentitySource.session:semantic,PressureActionIdIdentitySource.session_open:semantic,PressureActionIdIdentitySource.generation:semantic,PressureActionIdIdentitySource.client_key_digest:semantic",
    "source_bindings=PressureActionIdIdentitySource.governor_id>governor-id,PressureActionIdIdentitySource.session>session,PressureActionIdIdentitySource.session_open>session-open,PressureActionIdIdentitySource.generation>generation,PressureActionIdIdentitySource.client_key_digest>client-key-digest",
    "external_semantic_fields=identity-version,digest-domain,canonical-field-order",
    "semantic_fields=identity-version,digest-domain,canonical-field-order,governor-id,session,session-open,generation,client-key-digest",
    "excluded_fields=none",
    "consumers=Governor::apply_memory_pressure,PressureReceipt::action_id,DegradationEvent::pressure_action_id,buffered_degradation_event",
    "mutations=identity-version:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,digest-domain:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,canonical-field-order:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,governor-id:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,session:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,session-open:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,generation:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,client-key-digest:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_pressure_action_id_identity_fields",
    "transport_guard=session_identity_transport_is_current",
    "version_guard=crates/fs-session/src/governor.rs#governor_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-session:pressure-action-id",
];

/// Owner-local bounded submission-agent key declaration.
#[allow(dead_code)]
pub const SUBMISSION_AGENT_KEY_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-session:submission-agent-key",
    "version_const=SUBMISSION_AGENT_KEY_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-session.idempotency-agent.v1",
    "domain_const=IDEMPOTENCY_AGENT_DOMAIN",
    "encoder=submission_agent_key_identity",
    "encoder_helpers=none",
    "schema_constants=SUBMISSION_AGENT_KEY_IDENTITY_VERSION,IDEMPOTENCY_AGENT_DOMAIN,MAX_IDEMPOTENCY_INPUT_BYTES,crates/fs-session/src/governor/recovery.rs#TERMINAL_SCHEMA_VERSION",
    "schema_functions=bounded_request_digest,Governor::submission_request_id,Governor::idempotency_key,crates/fs-session/src/governor/recovery.rs#encode_submission_payload,crates/fs-session/src/governor/recovery.rs#decode_submission_payload,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=none",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=SubmissionAgentKeyIdentitySource",
    "source_fields=SubmissionAgentKeyIdentitySource.value:semantic",
    "source_bindings=SubmissionAgentKeyIdentitySource.value>agent-key-bytes",
    "external_semantic_fields=identity-version,digest-domain",
    "semantic_fields=identity-version,digest-domain,agent-key-bytes",
    "excluded_fields=none",
    "consumers=Governor::submission_request_id,Governor::idempotency_key,SubmissionRequestId::content_hash",
    "mutations=identity-version:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,digest-domain:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,agent-key-bytes:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_submission_agent_key_identity_fields",
    "transport_guard=session_identity_transport_is_current",
    "version_guard=crates/fs-session/src/governor.rs#governor_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-session:submission-agent-key",
];

/// Owner-local bounded canonical-program declaration.
#[allow(dead_code)]
pub const SUBMISSION_PROGRAM_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-session:submission-program",
    "version_const=SUBMISSION_PROGRAM_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-session.idempotency-program.v1",
    "domain_const=IDEMPOTENCY_PROGRAM_DOMAIN",
    "encoder=submission_program_identity",
    "encoder_helpers=none",
    "schema_constants=SUBMISSION_PROGRAM_IDENTITY_VERSION,IDEMPOTENCY_PROGRAM_DOMAIN,MAX_IDEMPOTENCY_INPUT_BYTES,crates/fs-session/src/governor/recovery.rs#TERMINAL_SCHEMA_VERSION",
    "schema_functions=bounded_request_digest,Governor::submission_request_id,Governor::idempotency_key,Governor::submit_once_durable,crates/fs-session/src/governor/recovery.rs#encode_submission_payload,crates/fs-session/src/governor/recovery.rs#decode_submission_payload,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=none",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=SubmissionProgramIdentitySource",
    "source_fields=SubmissionProgramIdentitySource.value:semantic",
    "source_bindings=SubmissionProgramIdentitySource.value>canonical-program-bytes",
    "external_semantic_fields=identity-version,digest-domain",
    "semantic_fields=identity-version,digest-domain,canonical-program-bytes",
    "excluded_fields=none",
    "consumers=Governor::submission_request_id,Governor::submit_once_durable,SubmissionRequestId",
    "mutations=identity-version:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,digest-domain:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,canonical-program-bytes:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_submission_program_identity_fields",
    "transport_guard=session_identity_transport_is_current",
    "version_guard=crates/fs-session/src/governor.rs#governor_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-session:submission-program",
];

/// Owner-local retryable submission authority declaration.
#[allow(dead_code)]
pub const SUBMISSION_REQUEST_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-session:submission-request-id",
    "version_const=SUBMISSION_REQUEST_IDENTITY_VERSION",
    "version=2",
    "domain=org.frankensim.fs-session.submission-request-id.v2",
    "domain_const=SUBMISSION_REQUEST_ID_DOMAIN",
    "encoder=submission_request_authority_identity",
    "encoder_helpers=identity_push_hash",
    "schema_constants=SUBMISSION_REQUEST_IDENTITY_VERSION,SUBMISSION_REQUEST_ID_DOMAIN,crates/fs-session/src/governor/recovery.rs#TERMINAL_SCHEMA_VERSION",
    "schema_functions=Governor::submission_request_id,Governor::submit_once,Governor::submit_once_durable,crates/fs-session/src/governor/recovery.rs#encode_submission_payload,crates/fs-session/src/governor/recovery.rs#decode_submission_payload,crates/fs-session/src/governor/recovery.rs#encode_submission_terminal_receipt,crates/fs-session/src/governor/recovery.rs#decode_submission_terminal_receipt,buffered_submission_success,buffered_submission_failure,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=fs-session:session-open-id,fs-session:submission-agent-key,fs-session:submission-program",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=SubmissionRequestIdIdentitySource",
    "source_fields=SubmissionRequestIdIdentitySource.governor_id:semantic,SubmissionRequestIdIdentitySource.session:semantic,SubmissionRequestIdIdentitySource.session_open:semantic,SubmissionRequestIdIdentitySource.generation:semantic,SubmissionRequestIdIdentitySource.key_hash:semantic,SubmissionRequestIdIdentitySource.request_hash:nonsemantic:checked-payload-not-slot-identity",
    "source_bindings=SubmissionRequestIdIdentitySource.governor_id>governor-id,SubmissionRequestIdIdentitySource.session>session,SubmissionRequestIdIdentitySource.session_open>session-open,SubmissionRequestIdIdentitySource.generation>generation,SubmissionRequestIdIdentitySource.key_hash>agent-key-hash",
    "external_semantic_fields=identity-version,digest-domain,canonical-field-order",
    "semantic_fields=identity-version,digest-domain,canonical-field-order,governor-id,session,session-open,generation,agent-key-hash",
    "excluded_fields=none",
    "consumers=Governor::submit_once,Governor::submit_once_durable,SubmissionReceipt,MeterReportId,fs-ledger::session_registry::SessionMutationClaim",
    "mutations=identity-version:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,digest-domain:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,canonical-field-order:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,governor-id:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,session:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,session-open:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,generation:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,agent-key-hash:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently",
    "nonsemantic_mutations=SubmissionRequestIdIdentitySource.request_hash:crates/fs-session/src/governor.rs#governor_identity_nonsemantic_fields_do_not_move",
    "field_guard=classify_submission_request_id_identity_fields",
    "transport_guard=session_identity_transport_is_current",
    "version_guard=crates/fs-session/src/governor.rs#governor_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-session:submission-request-id",
];

/// Owner-local retained-evidence declaration.
#[allow(dead_code)]
pub const RETAINED_EVIDENCE_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-session:retained-evidence",
    "version_const=RETAINED_EVIDENCE_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-session.retained-evidence.v1",
    "domain_const=RETAINED_EVIDENCE_DOMAIN",
    "encoder=retained_evidence_identity",
    "encoder_helpers=none",
    "schema_constants=RETAINED_EVIDENCE_IDENTITY_VERSION,RETAINED_EVIDENCE_DOMAIN,MAX_RETAINED_EVIDENCE_BYTES,MAX_CHECKPOINT_CLAIM_BYTES,crates/fs-session/src/governor/recovery.rs#TERMINAL_SCHEMA_VERSION",
    "schema_functions=RetainedEvidence::capture,utf8_prefix,crates/fs-session/src/governor/recovery.rs#encode_evidence,crates/fs-session/src/governor/recovery.rs#decode_evidence,evidence_json,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=none",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=RetainedEvidenceIdentitySource",
    "source_fields=RetainedEvidenceIdentitySource.value:semantic",
    "source_bindings=RetainedEvidenceIdentitySource.value>complete-evidence-bytes",
    "external_semantic_fields=identity-version,digest-domain",
    "semantic_fields=identity-version,digest-domain,complete-evidence-bytes",
    "excluded_fields=none",
    "consumers=SubmissionReceipt,DegradationEvent,PauseAcknowledgement,buffered_degradation_event",
    "mutations=identity-version:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,digest-domain:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,complete-evidence-bytes:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_retained_evidence_identity_fields",
    "transport_guard=session_identity_transport_is_current",
    "version_guard=crates/fs-session/src/governor.rs#governor_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-session:retained-evidence",
];

/// Owner-local pause-acknowledgement authority declaration.
#[allow(dead_code)]
pub const PAUSE_ACKNOWLEDGEMENT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-session:pause-acknowledgement-id",
    "version_const=PAUSE_ACKNOWLEDGEMENT_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-session.pause-acknowledgement-id.v1",
    "domain_const=PAUSE_ACKNOWLEDGEMENT_ID_DOMAIN",
    "encoder=pause_acknowledgement_authority_identity",
    "encoder_helpers=identity_push_hash",
    "schema_constants=PAUSE_ACKNOWLEDGEMENT_IDENTITY_VERSION,PAUSE_ACKNOWLEDGEMENT_ID_DOMAIN,crates/fs-session/src/governor/recovery.rs#TERMINAL_SCHEMA_VERSION",
    "schema_functions=crates/fs-session/src/governor/recovery.rs#pause_ack_authority,Governor::acknowledge_pause,crates/fs-session/src/governor/recovery.rs#encode_pause_ack_payload,crates/fs-session/src/governor/recovery.rs#decode_pause_ack_payload,crates/fs-session/src/governor/recovery.rs#encode_pause_ack_terminal_receipt,crates/fs-session/src/governor/recovery.rs#decode_pause_ack_terminal_receipt,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=fs-session:pressure-action-id",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=PauseAcknowledgementIdIdentitySource",
    "source_fields=PauseAcknowledgementIdIdentitySource.governor_id:semantic,PauseAcknowledgementIdIdentitySource.session:semantic,PauseAcknowledgementIdIdentitySource.gate_generation:semantic,PauseAcknowledgementIdIdentitySource.requested_ordinal:semantic",
    "source_bindings=PauseAcknowledgementIdIdentitySource.governor_id>governor-id,PauseAcknowledgementIdIdentitySource.session>session,PauseAcknowledgementIdIdentitySource.gate_generation>gate-generation,PauseAcknowledgementIdIdentitySource.requested_ordinal>requested-ordinal",
    "external_semantic_fields=identity-version,digest-domain,canonical-field-order",
    "semantic_fields=identity-version,digest-domain,canonical-field-order,governor-id,session,gate-generation,requested-ordinal",
    "excluded_fields=none",
    "consumers=Governor::acknowledge_pause,PauseAcknowledgement::request_id,DegradationEvent::pause_request_id,fs-ledger::session_registry::SessionMutationClaim",
    "mutations=identity-version:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,digest-domain:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,canonical-field-order:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,governor-id:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,session:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,gate-generation:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,requested-ordinal:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_pause_acknowledgement_id_identity_fields",
    "transport_guard=session_identity_transport_is_current",
    "version_guard=crates/fs-session/src/governor.rs#governor_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-session:pause-acknowledgement-id",
];

/// Owner-local resume-activation authority declaration.
#[allow(dead_code)]
pub const RESUME_ACTIVATION_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-session:resume-activation-id",
    "version_const=RESUME_ACTIVATION_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-session.resume-activation-id.v1",
    "domain_const=RESUME_ACTIVATION_ID_DOMAIN",
    "encoder=resume_activation_authority_identity",
    "encoder_helpers=identity_push_hash",
    "schema_constants=RESUME_ACTIVATION_IDENTITY_VERSION,RESUME_ACTIVATION_ID_DOMAIN,crates/fs-session/src/governor/recovery.rs#TERMINAL_SCHEMA_VERSION",
    "schema_functions=resume_activation_id,Governor::activate_resume,crates/fs-session/src/governor/recovery.rs#encode_activation_payload,crates/fs-session/src/governor/recovery.rs#decode_activation_payload,crates/fs-session/src/governor/recovery.rs#encode_activation_terminal_receipt,crates/fs-session/src/governor/recovery.rs#decode_activation_terminal_receipt,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=fs-session:gate-binding-id,fs-session:pause-acknowledgement-receipt,fs-session:session-open-id",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=ResumeActivationIdIdentitySource",
    "source_fields=ResumeActivationIdIdentitySource.governor_id:semantic,ResumeActivationIdIdentitySource.session:semantic,ResumeActivationIdIdentitySource.session_open:semantic,ResumeActivationIdIdentitySource.acknowledgement_hash:semantic,ResumeActivationIdIdentitySource.resume_generation:semantic",
    "source_bindings=ResumeActivationIdIdentitySource.governor_id>governor-id,ResumeActivationIdIdentitySource.session>session,ResumeActivationIdIdentitySource.session_open>session-open,ResumeActivationIdIdentitySource.acknowledgement_hash>acknowledgement-hash,ResumeActivationIdIdentitySource.resume_generation>resume-generation",
    "external_semantic_fields=identity-version,digest-domain,canonical-field-order",
    "semantic_fields=identity-version,digest-domain,canonical-field-order,governor-id,session,session-open,acknowledgement-hash,resume-generation",
    "excluded_fields=none",
    "consumers=Governor::activate_resume,ResumeActivationReceipt::activation_id,fs-ledger::session_registry::SessionMutationClaim",
    "mutations=identity-version:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,digest-domain:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,canonical-field-order:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,governor-id:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,session:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,session-open:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,acknowledgement-hash:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently,resume-generation:crates/fs-session/src/governor.rs#governor_authority_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_resume_activation_id_identity_fields",
    "transport_guard=session_identity_transport_is_current",
    "version_guard=crates/fs-session/src/governor.rs#governor_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-session:resume-activation-id",
];

/// Owner-local session-open receipt declaration.
#[allow(dead_code)]
pub const SESSION_OPEN_RECEIPT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-session:session-open-receipt",
    "version_const=SESSION_OPEN_RECEIPT_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-session.open-receipt.v1",
    "domain_const=SESSION_OPEN_RECEIPT_DOMAIN",
    "encoder=session_open_receipt_identity",
    "encoder_helpers=identity_push_framed,identity_push_hash",
    "schema_constants=SESSION_OPEN_RECEIPT_IDENTITY_VERSION,SESSION_OPEN_RECEIPT_DOMAIN,SESSION_OPEN_AUDIT_SCHEMA,crates/fs-session/src/governor/recovery.rs#TERMINAL_SCHEMA_VERSION",
    "schema_functions=session_open_receipt_hash,Governor::register_session,crates/fs-session/src/governor/recovery.rs#encode_open_payload,crates/fs-session/src/governor/recovery.rs#decode_open_payload,crates/fs-session/src/governor/recovery.rs#encode_open_receipt,crates/fs-session/src/governor/recovery.rs#decode_open_receipt,buffered_open_receipt,scoped_payload,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=fs-session:gate-binding-id,fs-session:session-open-id,fs-session:session-token-identity",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=SessionOpenReceiptIdentitySource",
    "source_fields=SessionOpenReceiptIdentitySource.open_id:semantic,SessionOpenReceiptIdentitySource.token_digest:semantic,SessionOpenReceiptIdentitySource.gate_identity:semantic,SessionOpenReceiptIdentitySource.ledger_scope:semantic",
    "source_bindings=SessionOpenReceiptIdentitySource.open_id>open-id,SessionOpenReceiptIdentitySource.token_digest>token-digest,SessionOpenReceiptIdentitySource.gate_identity>gate-presence+gate-identity,SessionOpenReceiptIdentitySource.ledger_scope>ledger-scope",
    "external_semantic_fields=identity-version,digest-domain,canonical-field-order,length-framing",
    "semantic_fields=identity-version,digest-domain,canonical-field-order,length-framing,open-id,token-digest,gate-presence,gate-identity,ledger-scope",
    "excluded_fields=none",
    "consumers=Governor::open_session,Governor::open_session_gated,ScopeFlushPermit,buffered_open_receipt,fs-ledger::session_registry::SessionTerminalRow",
    "mutations=identity-version:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,digest-domain:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,canonical-field-order:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,length-framing:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,open-id:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,token-digest:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,gate-presence:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,gate-identity:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,ledger-scope:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_session_open_receipt_identity_fields",
    "transport_guard=session_identity_transport_is_current",
    "version_guard=crates/fs-session/src/governor.rs#governor_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-session:session-open-receipt",
];

/// Owner-local meter receipt declaration.
#[allow(dead_code)]
pub const METER_RECEIPT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-session:meter-receipt",
    "version_const=METER_RECEIPT_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-session.meter-receipt.v1",
    "domain_const=METER_RECEIPT_DOMAIN",
    "encoder=meter_receipt_identity",
    "encoder_helpers=identity_push_enforcement,identity_push_framed,identity_push_hash,identity_push_snapshot",
    "schema_constants=METER_RECEIPT_IDENTITY_VERSION,METER_RECEIPT_DOMAIN,SESSION_METER_AUDIT_SCHEMA,crates/fs-session/src/governor/recovery.rs#TERMINAL_SCHEMA_VERSION",
    "schema_functions=meter_receipt_hash,push_meter_snapshot,push_enforcement_identity,Governor::commit_meter_locked,crates/fs-session/src/governor/recovery.rs#encode_meter_receipt,crates/fs-session/src/governor/recovery.rs#decode_meter_receipt,crates/fs-session/src/governor/recovery.rs#encode_meter_terminal_receipt,crates/fs-session/src/governor/recovery.rs#decode_meter_terminal_receipt,buffered_meter_receipt,enforcement_json,meter_snapshot_json,scoped_payload,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=fs-session:meter-report-id",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=MeterReceiptIdentitySource,ChargeIdentitySource,MeterSnapshotIdentitySource,EnforcementIdentitySource",
    "source_fields=MeterReceiptIdentitySource.report_id:semantic,MeterReceiptIdentitySource.commit_ordinal:semantic,MeterReceiptIdentitySource.delta:derived:nested-charge-fields-classified-separately,MeterReceiptIdentitySource.before:derived:nested-snapshot-fields-classified-separately,MeterReceiptIdentitySource.after:derived:nested-snapshot-fields-classified-separately,MeterReceiptIdentitySource.enforcement:derived:nested-enforcement-fields-classified-separately,ChargeIdentitySource.core_s_bits:semantic,ChargeIdentitySource.mem_peak_bytes:semantic,ChargeIdentitySource.wall_s_bits:semantic,MeterSnapshotIdentitySource.core_s_bits:semantic,MeterSnapshotIdentitySource.mem_peak_bytes:semantic,MeterSnapshotIdentitySource.wall_s_bits:semantic,MeterSnapshotIdentitySource.throttled:semantic,MeterSnapshotIdentitySource.paused:semantic,EnforcementIdentitySource.variant:semantic,EnforcementIdentitySource.resource:semantic,EnforcementIdentitySource.used_bits:semantic,EnforcementIdentitySource.granted_bits:semantic,EnforcementIdentitySource.resume_hint:semantic",
    "source_bindings=MeterReceiptIdentitySource.report_id>report-id,MeterReceiptIdentitySource.commit_ordinal>commit-ordinal,ChargeIdentitySource.core_s_bits>delta-core-seconds-bits,ChargeIdentitySource.mem_peak_bytes>delta-memory-bytes,ChargeIdentitySource.wall_s_bits>delta-wall-seconds-bits,MeterSnapshotIdentitySource.core_s_bits>before-core-seconds-bits+after-core-seconds-bits,MeterSnapshotIdentitySource.mem_peak_bytes>before-memory-bytes+after-memory-bytes,MeterSnapshotIdentitySource.wall_s_bits>before-wall-seconds-bits+after-wall-seconds-bits,MeterSnapshotIdentitySource.throttled>before-throttled+after-throttled,MeterSnapshotIdentitySource.paused>before-paused+after-paused,EnforcementIdentitySource.variant>enforcement-tag,EnforcementIdentitySource.resource>enforcement-resource,EnforcementIdentitySource.used_bits>enforcement-used-bits,EnforcementIdentitySource.granted_bits>enforcement-granted-bits,EnforcementIdentitySource.resume_hint>enforcement-resume-hint",
    "external_semantic_fields=identity-version,digest-domain,canonical-field-order",
    "semantic_fields=identity-version,digest-domain,canonical-field-order,report-id,commit-ordinal,delta-core-seconds-bits,delta-memory-bytes,delta-wall-seconds-bits,before-core-seconds-bits,after-core-seconds-bits,before-memory-bytes,after-memory-bytes,before-wall-seconds-bits,after-wall-seconds-bits,before-throttled,after-throttled,before-paused,after-paused,enforcement-tag,enforcement-resource,enforcement-used-bits,enforcement-granted-bits,enforcement-resume-hint",
    "excluded_fields=none",
    "consumers=Governor::charge,SubmitOutcome,SubmissionReceipt,buffered_meter_receipt,fs-ledger::session_registry::SessionTerminalRow",
    "mutations=identity-version:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,digest-domain:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,canonical-field-order:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,report-id:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,commit-ordinal:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,delta-core-seconds-bits:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,delta-memory-bytes:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,delta-wall-seconds-bits:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,before-core-seconds-bits:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,after-core-seconds-bits:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,before-memory-bytes:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,after-memory-bytes:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,before-wall-seconds-bits:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,after-wall-seconds-bits:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,before-throttled:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,after-throttled:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,before-paused:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,after-paused:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,enforcement-tag:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,enforcement-resource:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,enforcement-used-bits:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,enforcement-granted-bits:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,enforcement-resume-hint:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_meter_receipt_identity_fields",
    "transport_guard=session_identity_transport_is_current",
    "version_guard=crates/fs-session/src/governor.rs#governor_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-session:meter-receipt",
];

/// Owner-local pressure receipt declaration.
#[allow(dead_code)]
pub const PRESSURE_RECEIPT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-session:pressure-receipt",
    "version_const=PRESSURE_RECEIPT_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-session.pressure-receipt.v1",
    "domain_const=PRESSURE_RECEIPT_DOMAIN",
    "encoder=pressure_receipt_identity",
    "encoder_helpers=identity_push_framed,identity_push_hash,identity_push_pressure_event",
    "schema_constants=PRESSURE_RECEIPT_IDENTITY_VERSION,PRESSURE_RECEIPT_DOMAIN,SESSION_DEGRADATION_AUDIT_SCHEMA,crates/fs-session/src/governor/recovery.rs#TERMINAL_SCHEMA_VERSION",
    "schema_functions=pressure_receipt_hash,push_pressure_event_identity,Governor::apply_memory_pressure,crates/fs-session/src/governor/recovery.rs#encode_pressure_terminal_receipt,crates/fs-session/src/governor/recovery.rs#decode_pressure_terminal_receipt,crates/fs-session/src/governor/recovery.rs#encode_degradation_event,crates/fs-session/src/governor/recovery.rs#decode_degradation_event,buffered_degradation_event,scoped_payload,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=fs-session:pressure-action-id,fs-session:retained-evidence",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=PressureReceiptIdentitySource,PressureEventIdentitySource,RetainedEvidenceReceiptIdentitySource",
    "source_fields=PressureReceiptIdentitySource.action_id:semantic,PressureReceiptIdentitySource.level:semantic,PressureReceiptIdentitySource.events:semantic,PressureEventIdentitySource.session:semantic,PressureEventIdentitySource.step_tag:semantic,PressureEventIdentitySource.pressure_level:semantic,PressureEventIdentitySource.phase_tag:semantic,PressureEventIdentitySource.attribution:semantic,PressureEventIdentitySource.ordinal:semantic,PressureEventIdentitySource.requested_ordinal:semantic,PressureEventIdentitySource.checkpoint:derived:nested-evidence-fields-classified-separately,PressureEventIdentitySource.gate_generation:semantic,PressureEventIdentitySource.pause_request_id:nonsemantic:event-routing-authority-not-receipt-preimage,PressureEventIdentitySource.pressure_action_id:nonsemantic:event-routing-authority-already-bound-by-receipt,RetainedEvidenceReceiptIdentitySource.preview:nonsemantic:diagnostic-preview-not-identity,RetainedEvidenceReceiptIdentitySource.byte_len:semantic,RetainedEvidenceReceiptIdentitySource.digest:semantic",
    "source_bindings=PressureReceiptIdentitySource.action_id>action-id,PressureReceiptIdentitySource.level>level,PressureReceiptIdentitySource.events>event-count+event-order,PressureEventIdentitySource.session>event-session,PressureEventIdentitySource.step_tag>event-step-tag,PressureEventIdentitySource.pressure_level>event-pressure-level,PressureEventIdentitySource.phase_tag>event-phase-tag,PressureEventIdentitySource.attribution>event-attribution,PressureEventIdentitySource.ordinal>event-ordinal,PressureEventIdentitySource.requested_ordinal>requested-ordinal-presence+requested-ordinal,PressureEventIdentitySource.gate_generation>gate-generation-presence+gate-generation,RetainedEvidenceReceiptIdentitySource.byte_len>checkpoint-byte-length,RetainedEvidenceReceiptIdentitySource.digest>checkpoint-digest",
    "external_semantic_fields=identity-version,digest-domain,canonical-field-order,length-framing",
    "semantic_fields=identity-version,digest-domain,canonical-field-order,length-framing,action-id,level,event-count,event-order,event-session,event-step-tag,event-pressure-level,event-phase-tag,event-attribution,event-ordinal,requested-ordinal-presence,requested-ordinal,gate-generation-presence,gate-generation,checkpoint-byte-length,checkpoint-digest",
    "excluded_fields=none",
    "consumers=Governor::apply_memory_pressure,PressureReceipt::content_hash,buffered_degradation_event,fs-ledger::session_registry::SessionTerminalRow",
    "mutations=identity-version:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,digest-domain:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,canonical-field-order:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,length-framing:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,action-id:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,level:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,event-count:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,event-order:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,event-session:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,event-step-tag:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,event-pressure-level:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,event-phase-tag:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,event-attribution:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,event-ordinal:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,requested-ordinal-presence:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,requested-ordinal:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,gate-generation-presence:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,gate-generation:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,checkpoint-byte-length:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,checkpoint-digest:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently",
    "nonsemantic_mutations=PressureEventIdentitySource.pause_request_id:crates/fs-session/src/governor.rs#governor_identity_nonsemantic_fields_do_not_move,PressureEventIdentitySource.pressure_action_id:crates/fs-session/src/governor.rs#governor_identity_nonsemantic_fields_do_not_move,RetainedEvidenceReceiptIdentitySource.preview:crates/fs-session/src/governor.rs#governor_identity_nonsemantic_fields_do_not_move",
    "field_guard=classify_pressure_receipt_identity_fields",
    "transport_guard=session_identity_transport_is_current",
    "version_guard=crates/fs-session/src/governor.rs#governor_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-session:pressure-receipt",
];

/// Owner-local terminal submission receipt declaration.
#[allow(dead_code)]
pub const SUBMISSION_RECEIPT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-session:submission-receipt",
    "version_const=SUBMISSION_RECEIPT_IDENTITY_VERSION",
    "version=3",
    "domain=org.frankensim.fs-session.submission-receipt.v3",
    "domain_const=SUBMISSION_RECEIPT_DOMAIN",
    "encoder=submission_receipt_identity",
    "encoder_helpers=identity_push_framed,identity_push_hash",
    "schema_constants=SUBMISSION_RECEIPT_IDENTITY_VERSION,SUBMISSION_RECEIPT_DOMAIN,SESSION_SUBMISSION_AUDIT_SCHEMA,crates/fs-session/src/governor/recovery.rs#TERMINAL_SCHEMA_VERSION",
    "schema_functions=submission_receipt,SubmissionReceipt::matches_success,SubmissionReceipt::matches_failure,crates/fs-session/src/governor/recovery.rs#encode_submission_payload,crates/fs-session/src/governor/recovery.rs#decode_submission_payload,crates/fs-session/src/governor/recovery.rs#encode_submission_terminal_receipt,crates/fs-session/src/governor/recovery.rs#decode_submission_terminal_receipt,buffered_submission_success,buffered_submission_failure,scoped_payload,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=fs-session:meter-receipt,fs-session:retained-evidence,fs-session:submission-request-id",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=SubmissionReceiptIdentitySource,ChargeIdentitySource,RetainedEvidenceReceiptIdentitySource",
    "source_fields=SubmissionReceiptIdentitySource.variant:semantic,SubmissionReceiptIdentitySource.request_id:semantic,SubmissionReceiptIdentitySource.ledger_scope:semantic,SubmissionReceiptIdentitySource.admission_ordinal:semantic,SubmissionReceiptIdentitySource.charge:derived:nested-charge-fields-classified-separately,SubmissionReceiptIdentitySource.meter_receipt:semantic,SubmissionReceiptIdentitySource.evidence:derived:nested-evidence-fields-classified-separately,ChargeIdentitySource.core_s_bits:semantic,ChargeIdentitySource.mem_peak_bytes:semantic,ChargeIdentitySource.wall_s_bits:semantic,RetainedEvidenceReceiptIdentitySource.preview:nonsemantic:diagnostic-preview-not-identity,RetainedEvidenceReceiptIdentitySource.byte_len:semantic,RetainedEvidenceReceiptIdentitySource.digest:semantic",
    "source_bindings=SubmissionReceiptIdentitySource.variant>terminal-tag,SubmissionReceiptIdentitySource.request_id>request-id,SubmissionReceiptIdentitySource.ledger_scope>ledger-scope,SubmissionReceiptIdentitySource.admission_ordinal>admission-ordinal,SubmissionReceiptIdentitySource.meter_receipt>meter-receipt,ChargeIdentitySource.core_s_bits>charge-core-seconds-bits,ChargeIdentitySource.mem_peak_bytes>charge-memory-bytes,ChargeIdentitySource.wall_s_bits>charge-wall-seconds-bits,RetainedEvidenceReceiptIdentitySource.byte_len>evidence-byte-length,RetainedEvidenceReceiptIdentitySource.digest>evidence-digest",
    "external_semantic_fields=identity-version,digest-domain,canonical-field-order,length-framing",
    "semantic_fields=identity-version,digest-domain,canonical-field-order,length-framing,terminal-tag,request-id,ledger-scope,admission-ordinal,meter-receipt,charge-core-seconds-bits,charge-memory-bytes,charge-wall-seconds-bits,evidence-byte-length,evidence-digest",
    "excluded_fields=none",
    "consumers=SubmitOutcome,Governor::submit_once,Governor::submit_once_durable,buffered_submission_success,buffered_submission_failure,fs-ledger::session_registry::SessionTerminalRow",
    "mutations=identity-version:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,digest-domain:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,canonical-field-order:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,length-framing:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,terminal-tag:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,request-id:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,ledger-scope:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,admission-ordinal:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,meter-receipt:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,charge-core-seconds-bits:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,charge-memory-bytes:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,charge-wall-seconds-bits:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,evidence-byte-length:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,evidence-digest:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently",
    "nonsemantic_mutations=RetainedEvidenceReceiptIdentitySource.preview:crates/fs-session/src/governor.rs#governor_identity_nonsemantic_fields_do_not_move",
    "field_guard=classify_submission_receipt_identity_fields",
    "transport_guard=session_identity_transport_is_current",
    "version_guard=crates/fs-session/src/governor.rs#governor_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-session:submission-receipt",
];

/// Owner-local pause-acknowledgement receipt declaration.
#[allow(dead_code)]
pub const PAUSE_ACKNOWLEDGEMENT_RECEIPT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-session:pause-acknowledgement-receipt",
    "version_const=PAUSE_ACKNOWLEDGEMENT_RECEIPT_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-session.pause-acknowledgement-receipt.v1",
    "domain_const=PAUSE_ACKNOWLEDGEMENT_RECEIPT_DOMAIN",
    "encoder=pause_acknowledgement_receipt_identity",
    "encoder_helpers=identity_push_framed,identity_push_hash,identity_push_pressure_event",
    "schema_constants=PAUSE_ACKNOWLEDGEMENT_RECEIPT_IDENTITY_VERSION,PAUSE_ACKNOWLEDGEMENT_RECEIPT_DOMAIN,SESSION_DEGRADATION_AUDIT_SCHEMA,crates/fs-session/src/governor/recovery.rs#TERMINAL_SCHEMA_VERSION",
    "schema_functions=pause_acknowledgement_hash,Governor::acknowledge_pause,crates/fs-session/src/governor/recovery.rs#encode_pause_ack_payload,crates/fs-session/src/governor/recovery.rs#decode_pause_ack_payload,crates/fs-session/src/governor/recovery.rs#encode_pause_ack_terminal_receipt,crates/fs-session/src/governor/recovery.rs#decode_pause_ack_terminal_receipt,crates/fs-session/src/governor/recovery.rs#encode_degradation_event,crates/fs-session/src/governor/recovery.rs#decode_degradation_event,buffered_degradation_event,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=fs-session:gate-binding-id,fs-session:pause-acknowledgement-id,fs-session:retained-evidence",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=PauseAcknowledgementReceiptIdentitySource,PressureEventIdentitySource,RetainedEvidenceReceiptIdentitySource",
    "source_fields=PauseAcknowledgementReceiptIdentitySource.governor_id:semantic,PauseAcknowledgementReceiptIdentitySource.session:semantic,PauseAcknowledgementReceiptIdentitySource.gate_generation:semantic,PauseAcknowledgementReceiptIdentitySource.requested_ordinal:semantic,PauseAcknowledgementReceiptIdentitySource.event:derived:nested-event-fields-classified-separately,PauseAcknowledgementReceiptIdentitySource.resume_generation:semantic,PauseAcknowledgementReceiptIdentitySource.gate_binding:semantic,PressureEventIdentitySource.session:semantic,PressureEventIdentitySource.step_tag:semantic,PressureEventIdentitySource.pressure_level:semantic,PressureEventIdentitySource.phase_tag:semantic,PressureEventIdentitySource.attribution:semantic,PressureEventIdentitySource.ordinal:semantic,PressureEventIdentitySource.requested_ordinal:semantic,PressureEventIdentitySource.checkpoint:derived:nested-evidence-fields-classified-separately,PressureEventIdentitySource.gate_generation:semantic,PressureEventIdentitySource.pause_request_id:nonsemantic:event-routing-authority-not-receipt-preimage,PressureEventIdentitySource.pressure_action_id:nonsemantic:event-routing-authority-not-receipt-preimage,RetainedEvidenceReceiptIdentitySource.preview:nonsemantic:diagnostic-preview-not-identity,RetainedEvidenceReceiptIdentitySource.byte_len:semantic,RetainedEvidenceReceiptIdentitySource.digest:semantic",
    "source_bindings=PauseAcknowledgementReceiptIdentitySource.governor_id>governor-id,PauseAcknowledgementReceiptIdentitySource.session>request-session,PauseAcknowledgementReceiptIdentitySource.gate_generation>request-gate-generation,PauseAcknowledgementReceiptIdentitySource.requested_ordinal>request-ordinal,PauseAcknowledgementReceiptIdentitySource.resume_generation>resume-generation,PauseAcknowledgementReceiptIdentitySource.gate_binding>gate-binding,PressureEventIdentitySource.session>event-session,PressureEventIdentitySource.step_tag>event-step-tag,PressureEventIdentitySource.pressure_level>event-pressure-level,PressureEventIdentitySource.phase_tag>event-phase-tag,PressureEventIdentitySource.attribution>event-attribution,PressureEventIdentitySource.ordinal>event-ordinal,PressureEventIdentitySource.requested_ordinal>event-requested-ordinal-presence+event-requested-ordinal,PressureEventIdentitySource.gate_generation>event-gate-generation-presence+event-gate-generation,RetainedEvidenceReceiptIdentitySource.byte_len>checkpoint-byte-length,RetainedEvidenceReceiptIdentitySource.digest>checkpoint-digest",
    "external_semantic_fields=identity-version,digest-domain,canonical-field-order,length-framing",
    "semantic_fields=identity-version,digest-domain,canonical-field-order,length-framing,governor-id,request-session,request-gate-generation,request-ordinal,resume-generation,gate-binding,event-session,event-step-tag,event-pressure-level,event-phase-tag,event-attribution,event-ordinal,event-requested-ordinal-presence,event-requested-ordinal,event-gate-generation-presence,event-gate-generation,checkpoint-byte-length,checkpoint-digest",
    "excluded_fields=none",
    "consumers=Governor::acknowledge_pause,PauseAcknowledgement::content_hash,Governor::activate_resume,fs-ledger::session_registry::SessionTerminalRow",
    "mutations=identity-version:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,digest-domain:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,canonical-field-order:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,length-framing:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,governor-id:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,request-session:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,request-gate-generation:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,request-ordinal:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,resume-generation:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,gate-binding:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,event-session:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,event-step-tag:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,event-pressure-level:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,event-phase-tag:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,event-attribution:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,event-ordinal:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,event-requested-ordinal-presence:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,event-requested-ordinal:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,event-gate-generation-presence:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,event-gate-generation:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,checkpoint-byte-length:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,checkpoint-digest:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently",
    "nonsemantic_mutations=PressureEventIdentitySource.pause_request_id:crates/fs-session/src/governor.rs#governor_identity_nonsemantic_fields_do_not_move,PressureEventIdentitySource.pressure_action_id:crates/fs-session/src/governor.rs#governor_identity_nonsemantic_fields_do_not_move,RetainedEvidenceReceiptIdentitySource.preview:crates/fs-session/src/governor.rs#governor_identity_nonsemantic_fields_do_not_move",
    "field_guard=classify_pause_acknowledgement_receipt_identity_fields",
    "transport_guard=session_identity_transport_is_current",
    "version_guard=crates/fs-session/src/governor.rs#governor_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-session:pause-acknowledgement-receipt",
];

/// Owner-local resume-activation receipt declaration.
#[allow(dead_code)]
pub const RESUME_ACTIVATION_RECEIPT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-session:resume-activation-receipt",
    "version_const=RESUME_ACTIVATION_RECEIPT_IDENTITY_VERSION",
    "version=1",
    "domain=org.frankensim.fs-session.resume-activation-receipt.v1",
    "domain_const=RESUME_ACTIVATION_RECEIPT_DOMAIN",
    "encoder=resume_activation_receipt_identity",
    "encoder_helpers=identity_push_hash",
    "schema_constants=RESUME_ACTIVATION_RECEIPT_IDENTITY_VERSION,RESUME_ACTIVATION_RECEIPT_DOMAIN,crates/fs-session/src/governor/recovery.rs#TERMINAL_SCHEMA_VERSION",
    "schema_functions=resume_activation_receipt,Governor::activate_resume,crates/fs-session/src/governor/recovery.rs#encode_activation_payload,crates/fs-session/src/governor/recovery.rs#decode_activation_payload,crates/fs-session/src/governor/recovery.rs#encode_activation_terminal_receipt,crates/fs-session/src/governor/recovery.rs#decode_activation_terminal_receipt,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_dependencies=fs-session:gate-binding-id,fs-session:pause-acknowledgement-receipt,fs-session:resume-activation-id",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=ResumeActivationReceiptIdentitySource",
    "source_fields=ResumeActivationReceiptIdentitySource.activation_id:semantic,ResumeActivationReceiptIdentitySource.acknowledgement_hash:semantic,ResumeActivationReceiptIdentitySource.gate_binding:semantic",
    "source_bindings=ResumeActivationReceiptIdentitySource.activation_id>activation-id,ResumeActivationReceiptIdentitySource.acknowledgement_hash>acknowledgement-hash,ResumeActivationReceiptIdentitySource.gate_binding>gate-binding",
    "external_semantic_fields=identity-version,digest-domain,canonical-field-order",
    "semantic_fields=identity-version,digest-domain,canonical-field-order,activation-id,acknowledgement-hash,gate-binding",
    "excluded_fields=none",
    "consumers=Governor::activate_resume,ResumeActivationReceipt::content_hash,fs-ledger::session_registry::SessionTerminalRow",
    "mutations=identity-version:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,digest-domain:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,canonical-field-order:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,activation-id:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,acknowledgement-hash:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently,gate-binding:crates/fs-session/src/governor.rs#governor_receipt_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_resume_activation_receipt_identity_fields",
    "transport_guard=session_identity_transport_is_current",
    "version_guard=crates/fs-session/src/governor.rs#governor_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-session:resume-activation-receipt",
];
/// Maximum bytes hashed from either canonical-idempotency-key input.
pub const MAX_IDEMPOTENCY_INPUT_BYTES: usize = 1024 * 1024;
// Fixed-width conservative framing for t, three byte lengths, and the two
// optional-field discriminants in one persisted event row.
const FLUSH_ROW_FRAMING_BYTES: usize =
    core::mem::size_of::<i64>() + 3 * core::mem::size_of::<u64>() + 2 * core::mem::size_of::<u8>();
// Conservative claim + terminal + batch-member framing. The fs-ledger API
// performs the authoritative exact bound check before opening its transaction.
const FLUSH_TERMINAL_FRAMING_BYTES: usize = 512;

/// Maximum sessions admitted by one governor.
pub const MAX_SESSIONS_PER_GOVERNOR: usize = 4096;
/// Maximum sessions sharing one exact ledger scope.
pub const MAX_SESSIONS_PER_SCOPE: usize = 1024;
/// Maximum distinct idempotency keys retained for one session.
pub const MAX_IDEMPOTENCY_KEYS_PER_SESSION: usize = 4096;
/// Maximum distinct metering reports retained for one session.
pub const MAX_METER_REPORTS_PER_SESSION: usize = 8192;
/// Maximum distinct pressure actions retained for one session.
pub const MAX_PRESSURE_ACTIONS_PER_SESSION: usize = 4096;
/// Maximum degradation events retained in memory for one scope.
pub const MAX_DEGRADATION_EVENTS_PER_SCOPE: usize = 65_536;
/// Maximum UTF-8 bytes retained from caller-controlled diagnostic evidence.
pub const MAX_RETAINED_EVIDENCE_BYTES: usize = 16 * 1024;
/// Legacy upper bound retained for source compatibility.
///
/// Pause completion no longer accepts caller-supplied checkpoint claims; it
/// requires a fixed-width, ledger-verified `SolverCheckpointReceipt`.
pub const MAX_CHECKPOINT_CLAIM_BYTES: usize = 1024 * 1024;
/// Maximum caller-controlled payload bytes retained for one ledger scope.
pub const MAX_RETAINED_BYTES_PER_SCOPE: usize = 64 * 1024 * 1024;
/// Maximum caller-controlled payload bytes retained by one governor.
pub const MAX_RETAINED_BYTES_PER_GOVERNOR: usize = 256 * 1024 * 1024;
/// Maximum event rows emitted by one bounded flush call.
pub const MAX_FLUSH_ROWS: usize = 1024;
/// Maximum encoded event bytes emitted by one bounded flush call.
pub const MAX_FLUSH_ENCODED_BYTES: usize = 4 * 1024 * 1024;
/// Maximum degradation events returned by one page request.
pub const MAX_EVENT_PAGE_ROWS: usize = 1024;

const MAX_IDEMPOTENCY_TERMINAL_RETAINED_BYTES: usize = MAX_RETAINED_EVIDENCE_BYTES;
const MAX_PAUSE_COMPLETION_METADATA_BYTES: usize = 512;
const MAX_PAUSE_COMPLETION_RETAINED_BYTES: usize =
    MAX_RETAINED_EVIDENCE_BYTES + MAX_PAUSE_COMPLETION_METADATA_BYTES;
const MAX_METER_RECEIPT_RETAINED_BYTES: usize = 1024;
const PRESSURE_ACTION_RETAINED_BYTES: usize = 4 * core::mem::size_of::<u64>() + 64;
const OPEN_REQUEST_RETAINED_BYTES: usize = 4 * core::mem::size_of::<u64>() + 96;
const SUBMISSION_REQUEST_RETAINED_BYTES: usize = 3 * 32 + 64;

static NEXT_GOVERNOR_ID: AtomicU64 = AtomicU64::new(1);

fn ephemeral_governor_id() -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    payload.extend_from_slice(&std::process::id().to_le_bytes());
    payload.extend_from_slice(
        &std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos()
            .to_le_bytes(),
    );
    payload.extend_from_slice(
        &NEXT_GOVERNOR_ID
            .fetch_add(1, Ordering::Relaxed)
            .to_le_bytes(),
    );
    fs_blake3::hash_domain(EPHEMERAL_GOVERNOR_ID_DOMAIN, &payload)
}

fn utf8_prefix(value: &str, max_bytes: usize) -> String {
    let mut end = 0;
    for (index, ch) in value.char_indices() {
        let next = index + ch.len_utf8();
        if next > max_bytes {
            break;
        }
        end = next;
    }
    value[..end].to_string()
}

/// Bounded retained evidence for caller-controlled diagnostics.
///
/// The preview is UTF-8-safe and capped, while `byte_len` plus the
/// domain-separated digest bind the complete original input without retaining
/// it in governor state.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RetainedEvidence {
    preview: String,
    byte_len: usize,
    digest: fs_blake3::ContentHash,
}

impl RetainedEvidence {
    fn capture(value: &str) -> Self {
        Self {
            preview: utf8_prefix(value, MAX_RETAINED_EVIDENCE_BYTES),
            byte_len: value.len(),
            digest: fs_blake3::hash_domain(RETAINED_EVIDENCE_DOMAIN, value.as_bytes()),
        }
    }

    /// Bounded UTF-8-safe diagnostic prefix.
    #[must_use]
    pub fn preview(&self) -> &str {
        &self.preview
    }

    /// Exact byte length of the complete original evidence.
    #[must_use]
    pub const fn byte_len(&self) -> usize {
        self.byte_len
    }

    /// Digest of the complete original evidence.
    #[must_use]
    pub const fn digest(&self) -> fs_blake3::ContentHash {
        self.digest
    }
}

/// Explicit stable nonce used to reconstruct one governor identity against
/// the same physical Design Ledger after a process restart.
///
/// This is an identity input, not an authentication secret. The caller must
/// persist and version it alongside its session orchestration state.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DurableGovernorNonce(fs_blake3::ContentHash);

impl DurableGovernorNonce {
    /// Construct an explicit nonce from exactly 32 caller-persisted bytes.
    #[must_use]
    pub const fn from_bytes(bytes: [u8; 32]) -> Self {
        Self(fs_blake3::ContentHash(bytes))
    }

    /// Exact persisted nonce bytes.
    #[must_use]
    pub const fn as_bytes(self) -> [u8; 32] {
        self.0.0
    }
}

/// Unforgeable authority to flush one exact scope from one governor.
#[derive(Clone, PartialEq, Eq)]
pub struct ScopeFlushPermit {
    governor_id: fs_blake3::ContentHash,
    ledger_scope: String,
}

impl ScopeFlushPermit {
    /// Exact immutable scope carried by this permit.
    #[must_use]
    pub fn ledger_scope(&self) -> &str {
        &self.ledger_scope
    }
}

impl core::fmt::Debug for ScopeFlushPermit {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("ScopeFlushPermit")
            .field("ledger_scope", &self.ledger_scope)
            .finish_non_exhaustive()
    }
}

/// Opaque authority for one retryable session-open request.
///
/// The private fields bind the request to one governor and `SessionId`. The
/// caller supplies only a bounded request key; exact token and gate identity
/// are committed by the first successful open and checked on every replay.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SessionOpenId {
    governor_id: fs_blake3::ContentHash,
    session: SessionId,
    content_hash: fs_blake3::ContentHash,
}

impl SessionOpenId {
    /// Session this request is allowed to open.
    #[must_use]
    pub const fn session(self) -> SessionId {
        self.session
    }

    /// Domain-separated request identity.
    #[must_use]
    pub const fn content_hash(self) -> fs_blake3::ContentHash {
        self.content_hash
    }
}

/// Receipt for an admitted or exactly replayed session open.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SessionOpenReceipt {
    open_id: SessionOpenId,
    token_digest: fs_blake3::ContentHash,
    gate_identity: Option<fs_blake3::ContentHash>,
    permit: ScopeFlushPermit,
    content_hash: fs_blake3::ContentHash,
}

impl SessionOpenReceipt {
    /// Exact retry authority consumed by this open.
    #[must_use]
    pub const fn open_id(&self) -> SessionOpenId {
        self.open_id
    }

    /// Digest of the complete admitted capability token.
    #[must_use]
    pub const fn token_digest(&self) -> fs_blake3::ContentHash {
        self.token_digest
    }

    /// Governor-local identity of the bound cancellation gate, when gated.
    #[must_use]
    pub const fn gate_identity(&self) -> Option<fs_blake3::ContentHash> {
        self.gate_identity
    }

    /// Replayable permit for the exact immutable ledger scope.
    #[must_use]
    pub fn flush_permit(&self) -> ScopeFlushPermit {
        self.permit.clone()
    }

    /// Content identity binding authority, token, gate, and scope.
    #[must_use]
    pub const fn content_hash(&self) -> fs_blake3::ContentHash {
        self.content_hash
    }
}

/// Opaque authority for one retryable meter report.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MeterReportId {
    governor_id: fs_blake3::ContentHash,
    session: SessionId,
    session_open: fs_blake3::ContentHash,
    generation: u64,
    content_hash: fs_blake3::ContentHash,
}

impl MeterReportId {
    /// Session whose meter this report can mutate.
    #[must_use]
    pub const fn session(self) -> SessionId {
        self.session
    }

    /// Session execution generation captured when the authority was minted.
    #[must_use]
    pub const fn generation(self) -> u64 {
        self.generation
    }

    /// Domain-separated report identity.
    #[must_use]
    pub const fn content_hash(self) -> fs_blake3::ContentHash {
        self.content_hash
    }
}

/// Exact consumption state before or after one committed report.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeterSnapshot {
    /// Total core-seconds.
    pub core_s: f64,
    /// Peak resident bytes.
    pub mem_peak_bytes: u64,
    /// Total wall seconds.
    pub wall_s: f64,
    /// Number of committed throttling verdicts.
    pub throttled: u32,
    /// Number of committed pause verdicts.
    pub paused: u32,
}

/// Receipt for one atomically committed or exactly replayed meter report.
#[derive(Debug, Clone, PartialEq)]
pub struct MeterReceipt {
    report_id: MeterReportId,
    commit_ordinal: u64,
    delta: Charge,
    before: MeterSnapshot,
    after: MeterSnapshot,
    enforcement: Enforcement,
    content_hash: fs_blake3::ContentHash,
}

impl MeterReceipt {
    /// Report authority consumed by this commit.
    #[must_use]
    pub const fn report_id(&self) -> MeterReportId {
        self.report_id
    }

    /// Global causal meter-commit ordinal allocated with the charge.
    #[must_use]
    pub const fn commit_ordinal(&self) -> u64 {
        self.commit_ordinal
    }

    /// Exact-bit charge payload.
    #[must_use]
    pub const fn delta(&self) -> Charge {
        self.delta
    }

    /// Meter state immediately before the commit.
    #[must_use]
    pub const fn before(&self) -> MeterSnapshot {
        self.before
    }

    /// Meter state immediately after the commit.
    #[must_use]
    pub const fn after(&self) -> MeterSnapshot {
        self.after
    }

    /// Enforcement decision derived from `after`.
    #[must_use]
    pub const fn enforcement(&self) -> &Enforcement {
        &self.enforcement
    }

    /// Content identity binding authority, payload, causal order, and states.
    #[must_use]
    pub const fn content_hash(&self) -> fs_blake3::ContentHash {
        self.content_hash
    }
}

/// Opaque authority for one retryable declared pressure action.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PressureActionId {
    governor_id: fs_blake3::ContentHash,
    session: SessionId,
    session_open: fs_blake3::ContentHash,
    generation: u64,
    content_hash: fs_blake3::ContentHash,
}

impl PressureActionId {
    /// Session this action can mutate.
    #[must_use]
    pub const fn session(self) -> SessionId {
        self.session
    }

    /// Cancellation/execution generation captured at minting.
    #[must_use]
    pub const fn generation(self) -> u64 {
        self.generation
    }

    /// Domain-separated action identity.
    #[must_use]
    pub const fn content_hash(self) -> fs_blake3::ContentHash {
        self.content_hash
    }
}

/// Receipt for one committed or exactly replayed pressure action.
#[derive(Debug, Clone, PartialEq)]
pub struct PressureReceipt {
    action_id: PressureActionId,
    level: u8,
    events: Vec<DegradationEvent>,
    content_hash: fs_blake3::ContentHash,
}

impl PressureReceipt {
    /// Action authority consumed by this receipt.
    #[must_use]
    pub const fn action_id(&self) -> PressureActionId {
        self.action_id
    }

    /// Exact declared pressure level.
    #[must_use]
    pub const fn level(&self) -> u8 {
        self.level
    }

    /// Canonical event prefix committed exactly once by the action.
    #[must_use]
    pub fn events(&self) -> &[DegradationEvent] {
        &self.events
    }

    /// Content identity binding authority, level, and emitted events.
    #[must_use]
    pub const fn content_hash(&self) -> fs_blake3::ContentHash {
        self.content_hash
    }
}

/// Opaque authority for one admitted program submission.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SubmissionRequestId {
    governor_id: fs_blake3::ContentHash,
    session: SessionId,
    session_open: fs_blake3::ContentHash,
    generation: u64,
    key_hash: fs_blake3::ContentHash,
    request_hash: fs_blake3::ContentHash,
    content_hash: fs_blake3::ContentHash,
}

impl SubmissionRequestId {
    /// Session this request can submit into.
    #[must_use]
    pub const fn session(self) -> SessionId {
        self.session
    }

    /// Domain-separated identity of agent key plus canonical program.
    #[must_use]
    pub const fn content_hash(self) -> fs_blake3::ContentHash {
        self.content_hash
    }
}

/// Result of one bounded flush chunk.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FlushReport {
    /// Rows atomically appended by this call.
    pub appended_rows: usize,
    /// Typed terminal receipts atomically committed or verified by this call.
    pub committed_terminals: usize,
    /// Conservatively encoded bytes admitted to the batch.
    pub encoded_bytes: usize,
    /// More scoped state was dirty at return; call again with the same permit
    /// and ledger instance.
    pub remaining_dirty: bool,
}

/// One metering delta reported by the executor.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Charge {
    /// Core-seconds consumed.
    pub core_s: f64,
    /// Peak resident bytes observed during the interval.
    pub mem_peak_bytes: u64,
    /// Wall seconds elapsed.
    pub wall_s: f64,
}

/// The governor's enforcement verdict — always structured, never a kill.
#[derive(Debug, Clone, PartialEq)]
pub enum Enforcement {
    /// Within grants.
    Ok,
    /// At/over a grant: reduce concurrency; work continues.
    Throttled {
        /// Which grant bound (core-s / mem / wall).
        resource: &'static str,
        /// Consumed so far.
        used: f64,
        /// The grant.
        granted: f64,
    },
    /// Past the hard bound: checkpoint and stop; resumable by policy.
    Paused {
        /// Which grant bound.
        resource: &'static str,
        /// Consumed so far.
        used: f64,
        /// The grant.
        granted: f64,
        /// How to continue (teaching text).
        resume_hint: String,
    },
}

/// The declared degradation ladder — the ORDER is the contract.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DegradationStep {
    /// Spill the coldest arenas to disk.
    SpillColdArenas,
    /// Coarsen adaptive resolutions.
    CoarsenAdaptively,
    /// Checkpoint (SolverState) and stop; resume when pressure clears.
    PauseSerializeResume,
}

/// The ladder in its declared order.
pub const LADDER: [DegradationStep; 3] = [
    DegradationStep::SpillColdArenas,
    DegradationStep::CoarsenAdaptively,
    DegradationStep::PauseSerializeResume,
];

/// How far a ladder step has actually gotten (bead gp3.13): the ledger
/// distinguishes a synchronous action, a REQUEST awaiting the solver's
/// checkpoint, and the acknowledged completion — a pause that was never
/// acknowledged can never read as complete.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StepPhase {
    /// The governor declared an orchestration action that its owning subsystem
    /// must still perform and attest (spill/coarsen).
    Declared,
    /// Cancellation was requested on the session's OWN gate; the orchestrator
    /// has not yet acknowledged with a ledger checkpoint receipt.
    Requested,
    /// The owning orchestrator acknowledged a verified checkpoint receipt.
    Complete,
}

/// Opaque in-process authority for acknowledging one exact pause request.
///
/// Private fields bind the governor, session, gate generation, and request
/// ordinal so stale or cross-governor acknowledgements cannot complete a
/// different pause generation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PauseRequestId {
    governor_id: fs_blake3::ContentHash,
    session: SessionId,
    gate_generation: u64,
    requested_ordinal: i64,
}

impl PauseRequestId {
    /// Session whose pause this request controls.
    #[must_use]
    pub const fn session(self) -> SessionId {
        self.session
    }

    /// Cancellation-gate generation requested to drain.
    #[must_use]
    pub const fn gate_generation(self) -> u64 {
        self.gate_generation
    }

    /// Deterministic request-event ordinal.
    #[must_use]
    pub const fn requested_ordinal(self) -> i64 {
        self.requested_ordinal
    }

    /// Content-addressed authority that a solver checkpoint receipt must bind.
    ///
    /// This value is safe to persist or pass to `fs-ledger`; private request
    /// fields remain opaque so callers cannot forge a different pause request.
    #[must_use]
    pub fn checkpoint_authority(self) -> fs_blake3::ContentHash {
        recovery::pause_ack_authority(self)
    }
}

/// A ledgered degradation event.
#[derive(Debug, Clone, PartialEq)]
pub struct DegradationEvent {
    /// The affected session.
    pub session: SessionId,
    /// Which ladder step fired.
    pub step: DegradationStep,
    /// Pressure level (1..=3) that triggered it.
    pub pressure_level: u8,
    /// How far the step actually got (request vs acknowledged completion).
    pub phase: StepPhase,
    /// Attribution text (what was spilled/coarsened/paused).
    pub attribution: String,
    /// Logical event ordinal (deterministic; ledger `t`).
    pub ordinal: i64,
    /// Requested-event ordinal acknowledged by a completion event.
    pub requested_ordinal: Option<i64>,
    /// Bounded ledger checkpoint-receipt identity carried by completions.
    pub checkpoint: Option<RetainedEvidence>,
    /// Cancellation-gate generation for pause request/completion events.
    pub gate_generation: Option<u64>,
    /// Opaque acknowledgement authority for pause request/completion events.
    pub pause_request_id: Option<PauseRequestId>,
    /// Retry authority that committed this declared action.
    pub pressure_action_id: Option<PressureActionId>,
}

/// Successful pause acknowledgement plus the fresh gate that resumed work
/// must adopt and pass to [`Governor::activate_resume`]. The prior generation
/// remains permanently requested so old workers still drain; it is never reset
/// in place.
#[derive(Debug, Clone)]
pub struct PauseAcknowledgement {
    request_id: PauseRequestId,
    event: DegradationEvent,
    resume_gate: Arc<CancelGate>,
    resume_generation: u64,
    gate_binding: fs_blake3::ContentHash,
    content_hash: fs_blake3::ContentHash,
}

impl PauseAcknowledgement {
    /// Exact request consumed by this acknowledgement.
    #[must_use]
    pub const fn request_id(&self) -> PauseRequestId {
        self.request_id
    }

    /// Ledgerable completion event retained by the governor.
    #[must_use]
    pub const fn event(&self) -> &DegradationEvent {
        &self.event
    }

    /// Fresh unrequested gate for the next resumed generation.
    #[must_use]
    pub fn resume_gate(&self) -> Arc<CancelGate> {
        Arc::clone(&self.resume_gate)
    }

    /// Generation carried by the resume gate.
    #[must_use]
    pub const fn resume_generation(&self) -> u64 {
        self.resume_generation
    }

    /// Restart-stable semantic binding for the fresh gate generation. The
    /// process-local [`Arc`] may be rebound during recovery without changing
    /// this identity.
    #[must_use]
    pub const fn gate_binding(&self) -> fs_blake3::ContentHash {
        self.gate_binding
    }

    /// Content identity of the checkpoint acknowledgement terminal.
    #[must_use]
    pub const fn content_hash(&self) -> fs_blake3::ContentHash {
        self.content_hash
    }
}

/// Opaque authority for the idempotent activation of one acknowledged resume
/// generation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ResumeActivationId {
    governor_id: fs_blake3::ContentHash,
    session: SessionId,
    session_open: fs_blake3::ContentHash,
    resume_generation: u64,
    content_hash: fs_blake3::ContentHash,
}

impl ResumeActivationId {
    /// Session whose resumed generation this authority activates.
    #[must_use]
    pub const fn session(self) -> SessionId {
        self.session
    }

    /// Resumed gate generation named by this authority.
    #[must_use]
    pub const fn resume_generation(self) -> u64 {
        self.resume_generation
    }

    /// Domain-separated activation authority.
    #[must_use]
    pub const fn content_hash(self) -> fs_blake3::ContentHash {
        self.content_hash
    }
}

/// Structured receipt for an activated or exactly replayed resume generation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ResumeActivationReceipt {
    activation_id: ResumeActivationId,
    acknowledgement_hash: fs_blake3::ContentHash,
    gate_binding: fs_blake3::ContentHash,
    content_hash: fs_blake3::ContentHash,
}

impl ResumeActivationReceipt {
    /// Exact activation authority consumed by this receipt.
    #[must_use]
    pub const fn activation_id(self) -> ResumeActivationId {
        self.activation_id
    }

    /// Checkpoint acknowledgement this activation adopts.
    #[must_use]
    pub const fn acknowledgement_hash(self) -> fs_blake3::ContentHash {
        self.acknowledgement_hash
    }

    /// Restart-stable semantic gate binding adopted by workers.
    #[must_use]
    pub const fn gate_binding(self) -> fs_blake3::ContentHash {
        self.gate_binding
    }

    /// Content identity of the terminal activation receipt.
    #[must_use]
    pub const fn content_hash(self) -> fs_blake3::ContentHash {
        self.content_hash
    }
}

/// Opaque content identity for one terminal idempotent submission.
///
/// The private field prevents callers from minting receipts from arbitrary
/// integers. Identity binds the opaque request, admission order, immutable
/// ledger scope, and exact terminal meter receipt or failure diagnosis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SubmissionReceipt(fs_blake3::ContentHash);

impl SubmissionReceipt {
    /// Domain-separated content hash carried by this receipt.
    #[must_use]
    pub const fn content_hash(self) -> fs_blake3::ContentHash {
        self.0
    }

    /// Recompute and verify a successful terminal receipt.
    #[must_use]
    pub fn matches_success(
        self,
        request_id: SubmissionRequestId,
        ledger_scope: &str,
        admission_ordinal: u64,
        charge: Charge,
        meter_receipt: &MeterReceipt,
    ) -> bool {
        self == submission_receipt(
            request_id,
            ledger_scope,
            admission_ordinal,
            &SubmissionCompletion::Done(charge, meter_receipt.clone()),
        )
    }

    /// Recompute and verify a failed terminal receipt.
    #[must_use]
    pub fn matches_failure(
        self,
        request_id: SubmissionRequestId,
        ledger_scope: &str,
        admission_ordinal: u64,
        evidence: &RetainedEvidence,
    ) -> bool {
        self == submission_receipt(
            request_id,
            ledger_scope,
            admission_ordinal,
            &SubmissionCompletion::Failed(evidence.clone()),
        )
    }
}

impl core::fmt::Display for SubmissionReceipt {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        core::fmt::Display::fmt(&self.0, f)
    }
}

/// Outcome of an idempotent submission.
#[derive(Debug, Clone, PartialEq)]
pub enum SubmitOutcome {
    /// This call executed the work.
    Executed {
        /// Admission order reserved before caller work began.
        admission_ordinal: u64,
        /// The charge recorded.
        charge: Charge,
        /// Enforcement decision produced by committing that charge.
        enforcement: Enforcement,
        /// Causal pre/post meter receipt committed atomically with terminal publication.
        meter_receipt: MeterReceipt,
        /// Content-derived terminal receipt.
        receipt: SubmissionReceipt,
    },
    /// The key had already executed (or raced and lost): same receipt,
    /// NO additional charge.
    Duplicate {
        /// Original admission order.
        admission_ordinal: u64,
        /// The original execution's receipt.
        receipt: SubmissionReceipt,
        /// The original execution's enforcement decision.
        enforcement: Enforcement,
        /// The original execution's exact causal meter receipt.
        meter_receipt: MeterReceipt,
    },
    /// The one attempted execution failed before a charge could be committed.
    /// The key remains terminal: all duplicates receive this same receipt and
    /// diagnosis, and an explicit retry requires a new key.
    Failed {
        /// Original admission order.
        admission_ordinal: u64,
        /// The failed execution's receipt.
        receipt: SubmissionReceipt,
        /// Bounded preview plus full length/digest of the failure diagnosis.
        evidence: RetainedEvidence,
    },
    /// Another caller currently owns execution of this key. No waiting,
    /// execution, or charge occurred; poll/retry to observe its terminal state.
    InFlight,
    /// Rejected with guidance before execution.
    Refused(Box<Guidance>),
}

#[derive(Debug, Clone, Default)]
struct SessionMeters {
    core_s: f64,
    mem_peak_bytes: u64,
    wall_s: f64,
    throttled: u32,
    paused: u32,
}

impl SessionMeters {
    fn snapshot(&self) -> MeterSnapshot {
        MeterSnapshot {
            core_s: self.core_s,
            mem_peak_bytes: self.mem_peak_bytes,
            wall_s: self.wall_s,
            throttled: self.throttled,
            paused: self.paused,
        }
    }
}

#[derive(Debug, Clone)]
enum IdemState {
    Pending {
        admission_ordinal: u64,
        request_id: SubmissionRequestId,
        reserved_terminal_bytes: usize,
        reserved_meter_bytes: usize,
        durable_permit: Option<fs_ledger::session_registry::SessionClaimPermit>,
    },
    Done {
        admission_ordinal: u64,
        receipt: SubmissionReceipt,
        charge: Charge,
        meter_receipt: MeterReceipt,
        durable_permit: Option<fs_ledger::session_registry::SessionClaimPermit>,
    },
    Failed {
        admission_ordinal: u64,
        receipt: SubmissionReceipt,
        evidence: Arc<RetainedEvidence>,
        durable_permit: Option<fs_ledger::session_registry::SessionClaimPermit>,
    },
}

fn durable_submission_permit(
    state: &IdemState,
) -> Option<fs_ledger::session_registry::SessionClaimPermit> {
    match state {
        IdemState::Done { durable_permit, .. } | IdemState::Failed { durable_permit, .. } => {
            *durable_permit
        }
        IdemState::Pending { .. } => None,
    }
}

#[allow(clippy::large_enum_variant)] // Lock-local transition value; boxing would add a hot-path allocation.
enum SubmissionCompletion {
    Done(Charge, MeterReceipt),
    Failed(RetainedEvidence),
}

/// One meter-causal durable row. Successful submissions substitute their
/// self-contained terminal row for the private meter report that they own.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum DirtyCausalMutation {
    Meter(MeterReportId),
    Submission(SubmissionRequestId),
}

/// One indivisible durable control terminal. Variant order is causal when two
/// lifecycle terminals share the completion event ordinal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum DirtyControlMutation {
    Pressure(PressureActionId),
    PauseAcknowledgement(PauseRequestId),
    ResumeActivation(ResumeActivationId),
}

#[derive(Clone)]
struct OpenReplay {
    token_digest: fs_blake3::ContentHash,
    gate: Option<Arc<CancelGate>>,
    receipt: Arc<SessionOpenReceipt>,
}

#[derive(Clone, Copy)]
struct PressureReplay {
    level: u8,
    event_start: usize,
    event_len: usize,
    content_hash: fs_blake3::ContentHash,
}

struct BufferedLedgerEvent {
    session: [u8; 8],
    t: i64,
    kind: &'static str,
    payload: String,
}

impl BufferedLedgerEvent {
    fn as_row(&self) -> fs_ledger::EventRow<'_> {
        fs_ledger::EventRow {
            session: Some(self.session.as_slice()),
            t: self.t,
            kind: self.kind,
            payload: Some(&self.payload),
        }
    }

    fn encoded_len(&self) -> Result<usize, SessionError> {
        self.session
            .len()
            .checked_add(self.kind.len())
            .and_then(|bytes| bytes.checked_add(self.payload.len()))
            .and_then(|bytes| bytes.checked_add(FLUSH_ROW_FRAMING_BYTES))
            .ok_or(SessionError::LimitExceeded {
                resource: "flush_encoded_bytes",
                limit: MAX_FLUSH_ENCODED_BYTES,
                observed_at_least: usize::MAX,
            })
    }
}

fn buffered_open_receipt(
    ledger_scope: &str,
    open_id: SessionOpenId,
    receipt: &SessionOpenReceipt,
    token: &CapabilityToken,
) -> BufferedLedgerEvent {
    let gate_identity = receipt
        .gate_identity
        .map_or_else(|| "null".to_string(), |value| format!("\"{value}\""));
    BufferedLedgerEvent {
        session: open_id.session.0.to_be_bytes(),
        t: 0,
        kind: "session.open",
        payload: scoped_payload(
            SESSION_OPEN_AUDIT_SCHEMA,
            ledger_scope,
            &format!(
                "\"open_id\":\"{}\",\"token_digest\":\"{}\",\"gate_identity\":{gate_identity},\"receipt\":\"{}\",\"core_s_bits\":\"{:016x}\",\"mem_bytes\":{},\"wall_s_bits\":\"{:016x}\",\"cores\":{},\"ops\":{}",
                open_id.content_hash,
                receipt.token_digest,
                receipt.content_hash,
                token.core_s.to_bits(),
                token.mem_bytes,
                token.wall_s.to_bits(),
                token.cores,
                string_array_json(&token.ops),
            ),
        ),
    }
}

fn buffered_meter_receipt(
    ledger_scope: &str,
    report_id: MeterReportId,
    receipt: &MeterReceipt,
) -> Result<BufferedLedgerEvent, SessionError> {
    let before = receipt.before;
    let after = receipt.after;
    Ok(BufferedLedgerEvent {
        session: report_id.session.0.to_be_bytes(),
        t: i64::try_from(receipt.commit_ordinal).map_err(|_| SessionError::LimitExceeded {
            resource: "meter_commit_ordinal",
            limit: i64::MAX as usize,
            observed_at_least: usize::MAX,
        })?,
        kind: "session.meter-report",
        payload: scoped_payload(
            SESSION_METER_AUDIT_SCHEMA,
            ledger_scope,
            &format!(
                "\"session_open\":\"{}\",\"report_id\":\"{}\",\"generation\":{},\"commit_ordinal\":{},\"core_s_bits\":\"{:016x}\",\"mem_peak_bytes\":{},\"wall_s_bits\":\"{:016x}\",\"before\":{{\"core_s_bits\":\"{:016x}\",\"mem_peak_bytes\":{},\"wall_s_bits\":\"{:016x}\",\"throttled\":{},\"paused\":{}}},\"after\":{{\"core_s_bits\":\"{:016x}\",\"mem_peak_bytes\":{},\"wall_s_bits\":\"{:016x}\",\"throttled\":{},\"paused\":{}}},\"enforcement\":{},\"receipt\":\"{}\"",
                report_id.session_open,
                report_id.content_hash,
                report_id.generation,
                receipt.commit_ordinal,
                receipt.delta.core_s.to_bits(),
                receipt.delta.mem_peak_bytes,
                receipt.delta.wall_s.to_bits(),
                before.core_s.to_bits(),
                before.mem_peak_bytes,
                before.wall_s.to_bits(),
                before.throttled,
                before.paused,
                after.core_s.to_bits(),
                after.mem_peak_bytes,
                after.wall_s.to_bits(),
                after.throttled,
                after.paused,
                enforcement_json(&receipt.enforcement),
                receipt.content_hash,
            ),
        ),
    })
}

fn buffered_submission_success(
    ledger_scope: &str,
    request_id: SubmissionRequestId,
    state: &IdemState,
) -> Result<(BufferedLedgerEvent, (u64, SubmissionReceipt, u64)), SessionError> {
    let IdemState::Done {
        admission_ordinal,
        receipt,
        charge,
        meter_receipt,
        ..
    } = state
    else {
        return Err(SessionError::Persistence {
            what: format!(
                "causal submission index references non-success request {}",
                request_id.content_hash
            ),
        });
    };
    let derived_report_id = Governor::submission_meter_report_id(request_id);
    if meter_receipt.report_id != derived_report_id {
        return Err(SessionError::Persistence {
            what: format!(
                "submission {} terminal meter authority {} disagrees with derived authority {}",
                request_id.content_hash,
                meter_receipt.report_id.content_hash,
                derived_report_id.content_hash,
            ),
        });
    }
    let event_ordinal = meter_receipt.commit_ordinal;
    let session = request_id.session.0;
    let body = format!(
        "\"session\":{session},\"session_open\":\"{}\",\"generation\":{},\"request_id\":\"{}\",\"key_hash\":\"{}\",\"request_hash\":\"{}\",\"admission_ordinal\":{},\"meter_report_id\":\"{}\",\"meter_commit_ordinal\":{},\"meter_receipt\":\"{}\",\"receipt\":\"{receipt}\",\"core_s_bits\":\"{:016x}\",\"mem_peak_bytes\":{},\"wall_s_bits\":\"{:016x}\",\"before\":{},\"after\":{},\"enforcement\":{}",
        request_id.session_open,
        request_id.generation,
        request_id.content_hash,
        request_id.key_hash,
        request_id.request_hash,
        admission_ordinal,
        meter_receipt.report_id.content_hash,
        meter_receipt.commit_ordinal,
        meter_receipt.content_hash,
        charge.core_s.to_bits(),
        charge.mem_peak_bytes,
        charge.wall_s.to_bits(),
        meter_snapshot_json(meter_receipt.before),
        meter_snapshot_json(meter_receipt.after),
        enforcement_json(&meter_receipt.enforcement),
    );
    let event = BufferedLedgerEvent {
        session: session.to_be_bytes(),
        t: i64::try_from(event_ordinal).map_err(|_| SessionError::LimitExceeded {
            resource: "meter_commit_ordinal",
            limit: i64::MAX as usize,
            observed_at_least: usize::MAX,
        })?,
        kind: "session.idempotent-execution",
        payload: scoped_payload(SESSION_SUBMISSION_AUDIT_SCHEMA, ledger_scope, &body),
    };
    Ok((event, (*admission_ordinal, *receipt, event_ordinal)))
}

fn buffered_submission_failure(
    ledger_scope: &str,
    request_id: SubmissionRequestId,
    state: &IdemState,
) -> Result<(BufferedLedgerEvent, (u64, SubmissionReceipt, u64)), SessionError> {
    let IdemState::Failed {
        admission_ordinal,
        receipt,
        evidence,
        ..
    } = state
    else {
        return Err(SessionError::Persistence {
            what: format!(
                "failed-submission index references non-failure request {}",
                request_id.content_hash
            ),
        });
    };
    let session = request_id.session.0;
    let body = format!(
        "\"session\":{session},\"session_open\":\"{}\",\"generation\":{},\"request_id\":\"{}\",\"key_hash\":\"{}\",\"request_hash\":\"{}\",\"admission_ordinal\":{},\"receipt\":\"{receipt}\",\"error_evidence\":{}",
        request_id.session_open,
        request_id.generation,
        request_id.content_hash,
        request_id.key_hash,
        request_id.request_hash,
        admission_ordinal,
        evidence_json(evidence),
    );
    let event = BufferedLedgerEvent {
        session: session.to_be_bytes(),
        t: i64::try_from(*admission_ordinal).map_err(|_| SessionError::LimitExceeded {
            resource: "submission_ordinal",
            limit: i64::MAX as usize,
            observed_at_least: usize::MAX,
        })?,
        kind: "session.idempotent-failure",
        payload: scoped_payload(SESSION_SUBMISSION_AUDIT_SCHEMA, ledger_scope, &body),
    };
    Ok((event, (*admission_ordinal, *receipt, *admission_ordinal)))
}

fn buffered_degradation_event(
    ledger_scope: &str,
    event: &DegradationEvent,
    action_receipt_hash: fs_blake3::ContentHash,
) -> Result<BufferedLedgerEvent, SessionError> {
    let pressure_action_id = event
        .pressure_action_id
        .ok_or_else(|| SessionError::Persistence {
            what: format!(
                "degradation event {} lacks its pressure action authority",
                event.ordinal
            ),
        })?;
    let requested = event
        .requested_ordinal
        .map_or_else(|| "null".to_string(), |ordinal| ordinal.to_string());
    let checkpoint = event
        .checkpoint
        .as_ref()
        .map_or_else(|| "null".to_string(), evidence_json);
    let gate_generation = event
        .gate_generation
        .map_or_else(|| "null".to_string(), |value| value.to_string());
    Ok(BufferedLedgerEvent {
        session: event.session.0.to_be_bytes(),
        t: event.ordinal,
        kind: "session.degradation",
        payload: scoped_payload(
            SESSION_DEGRADATION_AUDIT_SCHEMA,
            ledger_scope,
            &format!(
                "\"session_open\":\"{}\",\"generation\":{},\"action_id\":\"{}\",\"action_receipt\":\"{}\",\"step\":\"{}\",\"level\":{},\"phase\":\"{}\",\"attribution\":\"{}\",\"requested_ordinal\":{requested},\"checkpoint\":{checkpoint},\"gate_generation\":{gate_generation}",
                pressure_action_id.session_open,
                pressure_action_id.generation,
                pressure_action_id.content_hash,
                action_receipt_hash,
                degradation_step_name(event.step),
                event.pressure_level,
                step_phase_name(event.phase),
                json_escape(&event.attribution),
            ),
        ),
    })
}

fn validate_resource(resource: &'static str, value: f64) -> Result<(), SessionError> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(SessionError::InvalidResource {
            resource,
            value,
            requirement: "must be finite and non-negative",
        })
    }
}

fn panic_evidence(payload: &(dyn std::any::Any + Send)) -> RetainedEvidence {
    if let Some(message) = payload.downcast_ref::<&str>() {
        RetainedEvidence::capture(message)
    } else if let Some(message) = payload.downcast_ref::<String>() {
        RetainedEvidence::capture(message)
    } else {
        RetainedEvidence::capture("submission work panicked with a non-string payload")
    }
}

fn push_framed(payload: &mut Vec<u8>, bytes: &[u8]) {
    payload.extend_from_slice(
        &u64::try_from(bytes.len())
            .expect("submission receipt field length fits u64")
            .to_le_bytes(),
    );
    payload.extend_from_slice(bytes);
}

fn bounded_request_digest(
    resource: &'static str,
    domain: &str,
    value: &str,
) -> Result<fs_blake3::ContentHash, SessionError> {
    if value.len() > MAX_IDEMPOTENCY_INPUT_BYTES {
        return Err(SessionError::LimitExceeded {
            resource,
            limit: MAX_IDEMPOTENCY_INPUT_BYTES,
            observed_at_least: value.len(),
        });
    }
    if value.trim().is_empty() {
        return Err(SessionError::Submission {
            what: format!("{resource} must be non-blank"),
        });
    }
    Ok(fs_blake3::hash_domain(domain, value.as_bytes()))
}

fn push_hash(payload: &mut Vec<u8>, hash: fs_blake3::ContentHash) {
    payload.extend_from_slice(hash.as_bytes());
}

fn capability_token_identity(token: &CapabilityToken) -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    payload.extend_from_slice(&token.session.0.to_le_bytes());
    payload.extend_from_slice(&token.core_s.to_bits().to_le_bytes());
    payload.extend_from_slice(&token.mem_bytes.to_le_bytes());
    payload.extend_from_slice(&token.wall_s.to_bits().to_le_bytes());
    payload.extend_from_slice(&token.cores.to_le_bytes());
    push_framed(&mut payload, token.ledger_scope.as_bytes());
    payload.extend_from_slice(
        &u64::try_from(token.ops.len())
            .expect("bounded operator count fits u64")
            .to_le_bytes(),
    );
    for operation in &token.ops {
        push_framed(&mut payload, operation.as_bytes());
    }
    fs_blake3::hash_domain(SESSION_TOKEN_IDENTITY_DOMAIN, &payload)
}

fn same_charge(left: Charge, right: Charge) -> bool {
    left.core_s.to_bits() == right.core_s.to_bits()
        && left.mem_peak_bytes == right.mem_peak_bytes
        && left.wall_s.to_bits() == right.wall_s.to_bits()
}

fn push_meter_snapshot(payload: &mut Vec<u8>, snapshot: MeterSnapshot) {
    payload.extend_from_slice(&snapshot.core_s.to_bits().to_le_bytes());
    payload.extend_from_slice(&snapshot.mem_peak_bytes.to_le_bytes());
    payload.extend_from_slice(&snapshot.wall_s.to_bits().to_le_bytes());
    payload.extend_from_slice(&snapshot.throttled.to_le_bytes());
    payload.extend_from_slice(&snapshot.paused.to_le_bytes());
}

fn meter_receipt_hash(
    report_id: MeterReportId,
    commit_ordinal: u64,
    delta: Charge,
    before: MeterSnapshot,
    after: MeterSnapshot,
    enforcement: &Enforcement,
) -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    push_hash(&mut payload, report_id.content_hash);
    payload.extend_from_slice(&commit_ordinal.to_le_bytes());
    payload.extend_from_slice(&delta.core_s.to_bits().to_le_bytes());
    payload.extend_from_slice(&delta.mem_peak_bytes.to_le_bytes());
    payload.extend_from_slice(&delta.wall_s.to_bits().to_le_bytes());
    push_meter_snapshot(&mut payload, before);
    push_meter_snapshot(&mut payload, after);
    push_enforcement_identity(&mut payload, enforcement);
    fs_blake3::hash_domain(METER_RECEIPT_DOMAIN, &payload)
}

fn push_pressure_event_identity(payload: &mut Vec<u8>, event: &DegradationEvent) {
    payload.extend_from_slice(&event.session.0.to_le_bytes());
    payload.push(match event.step {
        DegradationStep::SpillColdArenas => 0,
        DegradationStep::CoarsenAdaptively => 1,
        DegradationStep::PauseSerializeResume => 2,
    });
    payload.push(event.pressure_level);
    payload.push(match event.phase {
        StepPhase::Declared => 0,
        StepPhase::Requested => 1,
        StepPhase::Complete => 2,
    });
    push_framed(payload, event.attribution.as_bytes());
    payload.extend_from_slice(&event.ordinal.to_le_bytes());
    match event.requested_ordinal {
        Some(value) => {
            payload.push(1);
            payload.extend_from_slice(&value.to_le_bytes());
        }
        None => payload.push(0),
    }
    match &event.checkpoint {
        Some(checkpoint) => {
            payload.push(1);
            payload.extend_from_slice(
                &u64::try_from(checkpoint.byte_len)
                    .expect("bounded evidence fits u64")
                    .to_le_bytes(),
            );
            push_hash(payload, checkpoint.digest);
        }
        None => payload.push(0),
    }
    match event.gate_generation {
        Some(value) => {
            payload.push(1);
            payload.extend_from_slice(&value.to_le_bytes());
        }
        None => payload.push(0),
    }
}

fn pressure_receipt_hash(
    action_id: PressureActionId,
    level: u8,
    events: &[DegradationEvent],
) -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    push_hash(&mut payload, action_id.content_hash);
    payload.push(level);
    payload.extend_from_slice(
        &u64::try_from(events.len())
            .expect("bounded degradation event count fits u64")
            .to_le_bytes(),
    );
    for event in events {
        push_pressure_event_identity(&mut payload, event);
    }
    fs_blake3::hash_domain(PRESSURE_RECEIPT_DOMAIN, &payload)
}

fn session_open_receipt_hash(
    open_id: SessionOpenId,
    token_digest: fs_blake3::ContentHash,
    gate_identity: Option<fs_blake3::ContentHash>,
    ledger_scope: &str,
) -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    push_hash(&mut payload, open_id.content_hash);
    push_hash(&mut payload, token_digest);
    match gate_identity {
        Some(identity) => {
            payload.push(1);
            push_hash(&mut payload, identity);
        }
        None => payload.push(0),
    }
    push_framed(&mut payload, ledger_scope.as_bytes());
    fs_blake3::hash_domain(SESSION_OPEN_RECEIPT_DOMAIN, &payload)
}

fn session_gate_binding(open_id: SessionOpenId) -> fs_blake3::ContentHash {
    fs_blake3::hash_domain(GATE_BINDING_ID_DOMAIN, open_id.content_hash.as_bytes())
}

fn resumed_gate_binding(
    request_id: PauseRequestId,
    resume_generation: u64,
) -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    push_hash(&mut payload, request_id.governor_id);
    payload.extend_from_slice(&request_id.session.0.to_le_bytes());
    payload.extend_from_slice(&request_id.gate_generation.to_le_bytes());
    payload.extend_from_slice(&request_id.requested_ordinal.to_le_bytes());
    payload.extend_from_slice(&resume_generation.to_le_bytes());
    fs_blake3::hash_domain(GATE_BINDING_ID_DOMAIN, &payload)
}

fn pause_acknowledgement_hash(
    request_id: PauseRequestId,
    event: &DegradationEvent,
    resume_generation: u64,
    gate_binding: fs_blake3::ContentHash,
) -> fs_blake3::ContentHash {
    let mut payload = Vec::new();
    push_hash(&mut payload, request_id.governor_id);
    payload.extend_from_slice(&request_id.session.0.to_le_bytes());
    payload.extend_from_slice(&request_id.gate_generation.to_le_bytes());
    payload.extend_from_slice(&request_id.requested_ordinal.to_le_bytes());
    push_pressure_event_identity(&mut payload, event);
    payload.extend_from_slice(&resume_generation.to_le_bytes());
    push_hash(&mut payload, gate_binding);
    fs_blake3::hash_domain(PAUSE_ACKNOWLEDGEMENT_RECEIPT_DOMAIN, &payload)
}

fn resume_activation_id(
    governor_id: fs_blake3::ContentHash,
    session: SessionId,
    session_open: fs_blake3::ContentHash,
    acknowledgement_hash: fs_blake3::ContentHash,
    resume_generation: u64,
) -> ResumeActivationId {
    let mut payload = Vec::new();
    push_hash(&mut payload, governor_id);
    payload.extend_from_slice(&session.0.to_le_bytes());
    push_hash(&mut payload, session_open);
    push_hash(&mut payload, acknowledgement_hash);
    payload.extend_from_slice(&resume_generation.to_le_bytes());
    ResumeActivationId {
        governor_id,
        session,
        session_open,
        resume_generation,
        content_hash: fs_blake3::hash_domain(RESUME_ACTIVATION_ID_DOMAIN, &payload),
    }
}

fn resume_activation_receipt(
    activation_id: ResumeActivationId,
    acknowledgement_hash: fs_blake3::ContentHash,
    gate_binding: fs_blake3::ContentHash,
) -> ResumeActivationReceipt {
    let mut payload = Vec::new();
    push_hash(&mut payload, activation_id.content_hash);
    push_hash(&mut payload, acknowledgement_hash);
    push_hash(&mut payload, gate_binding);
    ResumeActivationReceipt {
        activation_id,
        acknowledgement_hash,
        gate_binding,
        content_hash: fs_blake3::hash_domain(RESUME_ACTIVATION_RECEIPT_DOMAIN, &payload),
    }
}

fn meter_transition(
    token: &CapabilityToken,
    before: &SessionMeters,
    delta: Charge,
) -> Result<(SessionMeters, Enforcement), SessionError> {
    validate_resource("core-seconds charge", delta.core_s)?;
    validate_resource("wall-seconds charge", delta.wall_s)?;
    let mut next = before.clone();
    let next_core_s = next.core_s + delta.core_s;
    let next_wall_s = next.wall_s + delta.wall_s;
    validate_resource("accumulated core-seconds", next_core_s)?;
    validate_resource("accumulated wall-seconds", next_wall_s)?;
    next.core_s = next_core_s;
    next.mem_peak_bytes = next.mem_peak_bytes.max(delta.mem_peak_bytes);
    next.wall_s = next_wall_s;
    let memory_past_hard = u128::from(next.mem_peak_bytes) * u128::from(HARD_FACTOR_DENOMINATOR)
        > u128::from(token.mem_bytes) * u128::from(HARD_FACTOR_NUMERATOR);
    #[allow(clippy::cast_precision_loss)]
    let memory_diagnostic = (next.mem_peak_bytes as f64, token.mem_bytes as f64);
    let hard_violation = if next.core_s > token.core_s * HARD_FACTOR {
        Some(("core-seconds", next.core_s, token.core_s))
    } else if memory_past_hard {
        Some(("memory-bytes", memory_diagnostic.0, memory_diagnostic.1))
    } else if next.wall_s > token.wall_s * HARD_FACTOR {
        Some(("wall-seconds", next.wall_s, token.wall_s))
    } else {
        None
    };
    let enforcement = if let Some((resource, used, granted)) = hard_violation {
        next.paused = next
            .paused
            .checked_add(1)
            .ok_or(SessionError::LimitExceeded {
                resource: "paused_meter_count",
                limit: u32::MAX as usize,
                observed_at_least: u32::MAX as usize,
            })?;
        Enforcement::Paused {
            resource,
            used,
            granted,
            resume_hint: format!(
                "checkpoint required before continuing; resume with a larger {resource} grant or \
                 a coarsened study — the caller must arrange and ledger the checkpoint explicitly"
            ),
        }
    } else {
        let throttle_violation = if next.core_s >= token.core_s {
            Some(("core-seconds", next.core_s, token.core_s))
        } else if next.mem_peak_bytes >= token.mem_bytes {
            Some(("memory-bytes", memory_diagnostic.0, memory_diagnostic.1))
        } else if next.wall_s >= token.wall_s {
            Some(("wall-seconds", next.wall_s, token.wall_s))
        } else {
            None
        };
        if let Some((resource, used, granted)) = throttle_violation {
            next.throttled = next
                .throttled
                .checked_add(1)
                .ok_or(SessionError::LimitExceeded {
                    resource: "throttled_meter_count",
                    limit: u32::MAX as usize,
                    observed_at_least: u32::MAX as usize,
                })?;
            Enforcement::Throttled {
                resource,
                used,
                granted,
            }
        } else {
            Enforcement::Ok
        }
    };
    Ok((next, enforcement))
}

fn push_enforcement_identity(payload: &mut Vec<u8>, enforcement: &Enforcement) {
    match enforcement {
        Enforcement::Ok => payload.push(0),
        Enforcement::Throttled {
            resource,
            used,
            granted,
        } => {
            payload.push(1);
            push_framed(payload, resource.as_bytes());
            payload.extend_from_slice(&used.to_bits().to_le_bytes());
            payload.extend_from_slice(&granted.to_bits().to_le_bytes());
        }
        Enforcement::Paused {
            resource,
            used,
            granted,
            resume_hint,
        } => {
            payload.push(2);
            push_framed(payload, resource.as_bytes());
            payload.extend_from_slice(&used.to_bits().to_le_bytes());
            payload.extend_from_slice(&granted.to_bits().to_le_bytes());
            push_framed(payload, resume_hint.as_bytes());
        }
    }
}

fn submission_receipt(
    request_id: SubmissionRequestId,
    ledger_scope: &str,
    admission_ordinal: u64,
    completion: &SubmissionCompletion,
) -> SubmissionReceipt {
    let mut payload = Vec::new();
    push_hash(&mut payload, request_id.content_hash);
    push_framed(&mut payload, ledger_scope.as_bytes());
    payload.extend_from_slice(&admission_ordinal.to_le_bytes());
    match completion {
        SubmissionCompletion::Done(charge, meter_receipt) => {
            payload.push(0);
            payload.extend_from_slice(&charge.core_s.to_bits().to_le_bytes());
            payload.extend_from_slice(&charge.mem_peak_bytes.to_le_bytes());
            payload.extend_from_slice(&charge.wall_s.to_bits().to_le_bytes());
            push_hash(&mut payload, meter_receipt.content_hash);
        }
        SubmissionCompletion::Failed(evidence) => {
            payload.push(1);
            payload.extend_from_slice(
                &u64::try_from(evidence.byte_len)
                    .expect("retained evidence length fits u64")
                    .to_le_bytes(),
            );
            payload.extend_from_slice(evidence.digest.as_bytes());
        }
    }
    SubmissionReceipt(fs_blake3::hash_domain(SUBMISSION_RECEIPT_DOMAIN, &payload))
}

fn evidence_json(evidence: &RetainedEvidence) -> String {
    format!(
        "{{\"preview\":\"{}\",\"byte_len\":{},\"digest\":\"{}\"}}",
        json_escape(&evidence.preview),
        evidence.byte_len,
        evidence.digest,
    )
}

fn json_escape(value: &str) -> String {
    use core::fmt::Write as _;

    let mut out = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c.is_control() => {
                let _ = write!(out, "\\u{:04x}", u32::from(c));
            }
            c => out.push(c),
        }
    }
    out
}

fn scoped_payload(schema: &str, ledger_scope: &str, body: &str) -> String {
    format!(
        "{{\"schema\":\"{}\",\"ledger_scope\":\"{}\",{body}}}",
        json_escape(schema),
        json_escape(ledger_scope),
    )
}

fn string_array_json(values: &[String]) -> String {
    let mut out = String::from("[");
    for (index, value) in values.iter().enumerate() {
        if index != 0 {
            out.push(',');
        }
        out.push('"');
        out.push_str(&json_escape(value));
        out.push('"');
    }
    out.push(']');
    out
}

fn enforcement_json(enforcement: &Enforcement) -> String {
    match enforcement {
        Enforcement::Ok => "{\"kind\":\"ok\"}".to_string(),
        Enforcement::Throttled {
            resource,
            used,
            granted,
        } => format!(
            "{{\"kind\":\"throttled\",\"resource\":\"{}\",\"used_bits\":\"{:016x}\",\"granted_bits\":\"{:016x}\"}}",
            json_escape(resource),
            used.to_bits(),
            granted.to_bits(),
        ),
        Enforcement::Paused {
            resource,
            used,
            granted,
            resume_hint,
        } => format!(
            "{{\"kind\":\"paused\",\"resource\":\"{}\",\"used_bits\":\"{:016x}\",\"granted_bits\":\"{:016x}\",\"resume_hint\":\"{}\"}}",
            json_escape(resource),
            used.to_bits(),
            granted.to_bits(),
            json_escape(resume_hint),
        ),
    }
}

fn meter_snapshot_json(snapshot: MeterSnapshot) -> String {
    format!(
        "{{\"core_s_bits\":\"{:016x}\",\"mem_peak_bytes\":{},\"wall_s_bits\":\"{:016x}\",\"throttled\":{},\"paused\":{}}}",
        snapshot.core_s.to_bits(),
        snapshot.mem_peak_bytes,
        snapshot.wall_s.to_bits(),
        snapshot.throttled,
        snapshot.paused,
    )
}

fn enforcement_retained_bytes(enforcement: &Enforcement) -> usize {
    match enforcement {
        Enforcement::Paused { resume_hint, .. } => resume_hint.len(),
        Enforcement::Ok | Enforcement::Throttled { .. } => 0,
    }
}

fn degradation_step_name(step: DegradationStep) -> &'static str {
    match step {
        DegradationStep::SpillColdArenas => "spill-cold-arenas",
        DegradationStep::CoarsenAdaptively => "coarsen-adaptively",
        DegradationStep::PauseSerializeResume => "pause-serialize-resume",
    }
}

fn step_phase_name(phase: StepPhase) -> &'static str {
    match phase {
        StepPhase::Declared => "declared",
        StepPhase::Requested => "requested",
        StepPhase::Complete => "complete",
    }
}

fn degradation_attribution(step: DegradationStep) -> &'static str {
    match step {
        DegradationStep::SpillColdArenas => {
            "declared spill of coldest arenas (least-recently-touched first)"
        }
        DegradationStep::CoarsenAdaptively => {
            "declared adaptive coarsening outside protected bands"
        }
        DegradationStep::PauseSerializeResume => {
            "requested pause on the session-owned gate: owning orchestrator must drain at a tile \
             boundary, persist a solver-state artifact, and obtain executor plus ledger receipts \
             before acknowledge_pause"
        }
    }
}

fn degradation_event_retained_bytes(event: &DegradationEvent) -> Result<usize, SessionError> {
    event
        .attribution
        .len()
        .checked_add(
            event
                .checkpoint
                .as_ref()
                .map_or(0, |checkpoint| checkpoint.preview.len()),
        )
        .ok_or(SessionError::LimitExceeded {
            resource: "retained_bytes_per_scope",
            limit: MAX_RETAINED_BYTES_PER_SCOPE,
            observed_at_least: usize::MAX,
        })
}

struct PreparedFlush {
    reservation_id: u64,
    generation: i64,
    revision: u64,
    next_flush_lane: u8,
    terminals: Vec<BufferedTerminal>,
    encoded_bytes: usize,
    open_marks: Vec<(SessionOpenId, fs_blake3::ContentHash)>,
    meter_report_marks: Vec<(MeterReportId, fs_blake3::ContentHash)>,
    idempotency_marks: Vec<(SubmissionRequestId, (u64, SubmissionReceipt, u64))>,
    control_marks: Vec<(i64, DirtyControlMutation, fs_blake3::ContentHash, usize)>,
}

struct FlushSnapshot {
    reservation_id: u64,
    generation: i64,
    revision: u64,
    next_flush_lane: u8,
    sources: Vec<FlushSource>,
}

#[allow(clippy::large_enum_variant)] // Boxing would add up to 1,024 allocations under the hot-path lock.
enum FlushSource {
    Open {
        open_id: SessionOpenId,
        receipt: Arc<SessionOpenReceipt>,
        token: Arc<CapabilityToken>,
    },
    Meter {
        indexed_ordinal: u64,
        report_id: MeterReportId,
        receipt: MeterReceipt,
    },
    Submission {
        indexed_ordinal: u64,
        request_id: SubmissionRequestId,
        state: IdemState,
        failed_lane: bool,
    },
    Pressure {
        indexed_ordinal: i64,
        action_id: PressureActionId,
        replay: PressureReplay,
        events: Vec<Arc<DegradationEvent>>,
    },
    PauseAcknowledgement {
        indexed_ordinal: i64,
        request_id: PauseRequestId,
        replay: PauseAcknowledgementReplay,
        event: Arc<DegradationEvent>,
        action_hash: fs_blake3::ContentHash,
        gate: Arc<CancelGate>,
        session_open: fs_blake3::ContentHash,
    },
    ResumeActivation {
        indexed_ordinal: i64,
        activation_id: ResumeActivationId,
        receipt: ResumeActivationReceipt,
    },
}

enum FlushMark {
    Open(SessionOpenId, fs_blake3::ContentHash),
    Meter(MeterReportId, fs_blake3::ContentHash),
    Submission(SubmissionRequestId, (u64, SubmissionReceipt, u64)),
    Control(i64, DirtyControlMutation, fs_blake3::ContentHash, usize),
}

struct BufferedTerminal {
    authority: fs_blake3::ContentHash,
    session_open: fs_blake3::ContentHash,
    kind: &'static str,
    session: SessionId,
    generation: u64,
    causal_ordinal: Option<u64>,
    payload: Vec<u8>,
    receipt: Vec<u8>,
    events: Vec<BufferedLedgerEvent>,
    permit: Option<fs_ledger::session_registry::SessionClaimPermit>,
}

impl BufferedTerminal {
    fn encoded_len(&self, ledger_scope: &str) -> Result<usize, SessionError> {
        let mut bytes = FLUSH_TERMINAL_FRAMING_BYTES
            .checked_add(self.kind.len())
            .and_then(|value| value.checked_add(ledger_scope.len()))
            .and_then(|value| value.checked_add(self.payload.len()))
            .and_then(|value| value.checked_add(self.receipt.len()))
            .ok_or(SessionError::LimitExceeded {
                resource: "flush_encoded_bytes",
                limit: MAX_FLUSH_ENCODED_BYTES,
                observed_at_least: usize::MAX,
            })?;
        for event in &self.events {
            bytes = bytes
                .checked_add(event.encoded_len()?)
                .ok_or(SessionError::LimitExceeded {
                    resource: "flush_encoded_bytes",
                    limit: MAX_FLUSH_ENCODED_BYTES,
                    observed_at_least: usize::MAX,
                })?;
        }
        Ok(bytes)
    }
}

fn push_bounded_terminal(
    terminals: &mut Vec<BufferedTerminal>,
    event_rows: &mut usize,
    encoded_bytes: &mut usize,
    ledger_scope: &str,
    terminal: BufferedTerminal,
) -> Result<bool, SessionError> {
    let terminal_bytes = terminal.encoded_len(ledger_scope)?;
    let next_events =
        event_rows
            .checked_add(terminal.events.len())
            .ok_or(SessionError::LimitExceeded {
                resource: "flush_rows",
                limit: MAX_FLUSH_ROWS,
                observed_at_least: usize::MAX,
            })?;
    let next_bytes =
        encoded_bytes
            .checked_add(terminal_bytes)
            .ok_or(SessionError::LimitExceeded {
                resource: "flush_encoded_bytes",
                limit: MAX_FLUSH_ENCODED_BYTES,
                observed_at_least: usize::MAX,
            })?;
    if terminal.events.len() > MAX_FLUSH_ROWS || terminal_bytes > MAX_FLUSH_ENCODED_BYTES {
        return Err(SessionError::LimitExceeded {
            resource: "flush_terminal_encoded_bytes",
            limit: MAX_FLUSH_ENCODED_BYTES,
            observed_at_least: terminal_bytes,
        });
    }
    if terminals.len() == MAX_FLUSH_ROWS
        || next_events > MAX_FLUSH_ROWS
        || next_bytes > MAX_FLUSH_ENCODED_BYTES
    {
        return Ok(false);
    }
    terminals.push(terminal);
    *event_rows = next_events;
    *encoded_bytes = next_bytes;
    Ok(true)
}

impl FlushSnapshot {
    #[allow(clippy::too_many_lines)] // One explicit off-lock materialization grammar for every terminal kind.
    fn materialize(self, ledger_scope: &str) -> Result<PreparedFlush, SessionError> {
        let mut terminals = Vec::with_capacity(self.sources.len().min(64));
        let mut event_rows = 0usize;
        let mut encoded_bytes = 0usize;
        let mut open_marks = Vec::new();
        let mut meter_report_marks = Vec::new();
        let mut idempotency_marks = Vec::new();
        let mut control_marks = Vec::new();

        for source in self.sources {
            let (terminal, mark) = match source {
                FlushSource::Open {
                    open_id,
                    receipt,
                    token,
                } => {
                    if capability_token_identity(&token) != receipt.token_digest {
                        return Err(SessionError::Persistence {
                            what: format!(
                                "session {} token no longer matches its immutable open receipt",
                                open_id.session.0
                            ),
                        });
                    }
                    let row = buffered_open_receipt(ledger_scope, open_id, &receipt, &token);
                    (
                        BufferedTerminal {
                            authority: open_id.content_hash,
                            session_open: receipt.content_hash,
                            kind: recovery::KIND_OPEN,
                            session: open_id.session,
                            generation: 0,
                            causal_ordinal: None,
                            payload: recovery::encode_open_payload(&token, receipt.gate_identity),
                            receipt: recovery::encode_open_receipt(&receipt),
                            events: vec![row],
                            permit: None,
                        },
                        FlushMark::Open(open_id, receipt.content_hash),
                    )
                }
                FlushSource::Meter {
                    indexed_ordinal,
                    report_id,
                    receipt,
                } => {
                    if receipt.commit_ordinal != indexed_ordinal {
                        return Err(SessionError::Persistence {
                            what: format!(
                                "meter report {} causal index ordinal {indexed_ordinal} disagrees with receipt ordinal {}",
                                report_id.content_hash, receipt.commit_ordinal
                            ),
                        });
                    }
                    let row = buffered_meter_receipt(ledger_scope, report_id, &receipt)?;
                    let content_hash = receipt.content_hash;
                    (
                        BufferedTerminal {
                            authority: report_id.content_hash,
                            session_open: report_id.session_open,
                            kind: recovery::KIND_METER,
                            session: report_id.session,
                            generation: report_id.generation,
                            causal_ordinal: Some(receipt.commit_ordinal),
                            payload: recovery::encode_meter_payload(receipt.delta),
                            receipt: recovery::encode_meter_terminal_receipt(&receipt),
                            events: vec![row],
                            permit: None,
                        },
                        FlushMark::Meter(report_id, content_hash),
                    )
                }
                FlushSource::Submission {
                    indexed_ordinal,
                    request_id,
                    state,
                    failed_lane,
                } => {
                    let (row, generation) = if failed_lane {
                        buffered_submission_failure(ledger_scope, request_id, &state)?
                    } else {
                        buffered_submission_success(ledger_scope, request_id, &state)?
                    };
                    if generation.2 != indexed_ordinal {
                        return Err(SessionError::Persistence {
                            what: format!(
                                "submission {} index ordinal {indexed_ordinal} disagrees with terminal ordinal {}",
                                request_id.content_hash, generation.2
                            ),
                        });
                    }
                    (
                        BufferedTerminal {
                            authority: request_id.content_hash,
                            session_open: request_id.session_open,
                            kind: recovery::KIND_SUBMISSION,
                            session: request_id.session,
                            generation: request_id.generation,
                            causal_ordinal: Some(generation.0),
                            payload: recovery::encode_submission_payload(request_id),
                            receipt: recovery::encode_submission_terminal_receipt(&state)?,
                            events: vec![row],
                            permit: durable_submission_permit(&state),
                        },
                        FlushMark::Submission(request_id, generation),
                    )
                }
                FlushSource::Pressure {
                    indexed_ordinal,
                    action_id,
                    replay,
                    events,
                } => {
                    if events.first().map(|event| event.ordinal) != Some(indexed_ordinal) {
                        return Err(SessionError::Persistence {
                            what: format!(
                                "pressure action {} control ordinal disagrees with its first event",
                                action_id.content_hash
                            ),
                        });
                    }
                    let owned_events: Vec<_> =
                        events.iter().map(|event| event.as_ref().clone()).collect();
                    let receipt = PressureReceipt {
                        action_id,
                        level: replay.level,
                        events: owned_events,
                        content_hash: replay.content_hash,
                    };
                    let rows = events
                        .iter()
                        .map(|event| {
                            buffered_degradation_event(ledger_scope, event, replay.content_hash)
                        })
                        .collect::<Result<Vec<_>, _>>()?;
                    let owned_event_count = events.len();
                    (
                        BufferedTerminal {
                            authority: action_id.content_hash,
                            session_open: action_id.session_open,
                            kind: recovery::KIND_PRESSURE,
                            session: action_id.session,
                            generation: action_id.generation,
                            causal_ordinal: Some(u64::try_from(indexed_ordinal).map_err(|_| {
                                SessionError::Persistence {
                                    what: format!(
                                        "control terminal has negative causal ordinal {indexed_ordinal}"
                                    ),
                                }
                            })?),
                            payload: recovery::encode_pressure_payload(replay.level),
                            receipt: recovery::encode_pressure_terminal_receipt(&receipt),
                            events: rows,
                            permit: None,
                        },
                        FlushMark::Control(
                            indexed_ordinal,
                            DirtyControlMutation::Pressure(action_id),
                            replay.content_hash,
                            owned_event_count,
                        ),
                    )
                }
                FlushSource::PauseAcknowledgement {
                    indexed_ordinal,
                    request_id,
                    replay,
                    event,
                    action_hash,
                    gate,
                    session_open,
                } => {
                    if event.ordinal != indexed_ordinal {
                        return Err(SessionError::Persistence {
                            what: "pause acknowledgement control ordinal disagrees with its completion event"
                                .to_string(),
                        });
                    }
                    let acknowledgement = PauseAcknowledgement {
                        request_id,
                        event: event.as_ref().clone(),
                        resume_gate: gate,
                        resume_generation: replay.resume_generation,
                        gate_binding: replay.gate_binding,
                        content_hash: replay.content_hash,
                    };
                    let evidence =
                        event
                            .checkpoint
                            .as_ref()
                            .ok_or_else(|| SessionError::Persistence {
                                what: "pause acknowledgement lacks checkpoint evidence".to_string(),
                            })?;
                    let row = buffered_degradation_event(ledger_scope, &event, action_hash)?;
                    (
                        BufferedTerminal {
                            authority: recovery::pause_ack_authority(request_id),
                            session_open,
                            kind: recovery::KIND_PAUSE_ACK,
                            session: request_id.session,
                            generation: replay.resume_generation,
                            causal_ordinal: Some(u64::try_from(indexed_ordinal).map_err(|_| {
                                SessionError::Persistence {
                                    what: format!(
                                        "control terminal has negative causal ordinal {indexed_ordinal}"
                                    ),
                                }
                            })?),
                            payload: recovery::encode_pause_ack_payload(evidence),
                            receipt: recovery::encode_pause_ack_terminal_receipt(
                                &acknowledgement,
                            ),
                            events: vec![row],
                            permit: None,
                        },
                        FlushMark::Control(
                            indexed_ordinal,
                            DirtyControlMutation::PauseAcknowledgement(request_id),
                            replay.content_hash,
                            1,
                        ),
                    )
                }
                FlushSource::ResumeActivation {
                    indexed_ordinal,
                    activation_id,
                    receipt,
                } => (
                    BufferedTerminal {
                        authority: activation_id.content_hash,
                        session_open: activation_id.session_open,
                        kind: recovery::KIND_RESUME_ACTIVATION,
                        session: activation_id.session,
                        generation: activation_id.resume_generation,
                        causal_ordinal: Some(u64::try_from(indexed_ordinal).map_err(|_| {
                            SessionError::Persistence {
                                what: format!(
                                    "control terminal has negative causal ordinal {indexed_ordinal}"
                                ),
                            }
                        })?),
                        payload: recovery::encode_activation_payload_parts(
                            receipt.acknowledgement_hash,
                            activation_id.resume_generation,
                            receipt.gate_binding,
                        ),
                        receipt: recovery::encode_activation_terminal_receipt(receipt),
                        events: Vec::new(),
                        permit: None,
                    },
                    FlushMark::Control(
                        indexed_ordinal,
                        DirtyControlMutation::ResumeActivation(activation_id),
                        receipt.content_hash,
                        0,
                    ),
                ),
            };

            if !push_bounded_terminal(
                &mut terminals,
                &mut event_rows,
                &mut encoded_bytes,
                ledger_scope,
                terminal,
            )? {
                break;
            }
            match mark {
                FlushMark::Open(open_id, content_hash) => {
                    open_marks.push((open_id, content_hash));
                }
                FlushMark::Meter(report_id, content_hash) => {
                    meter_report_marks.push((report_id, content_hash));
                }
                FlushMark::Submission(request_id, generation) => {
                    idempotency_marks.push((request_id, generation));
                }
                FlushMark::Control(ordinal, mutation, content_hash, owned_events) => {
                    control_marks.push((ordinal, mutation, content_hash, owned_events));
                }
            }
        }

        Ok(PreparedFlush {
            reservation_id: self.reservation_id,
            generation: self.generation,
            revision: self.revision,
            next_flush_lane: self.next_flush_lane,
            terminals,
            encoded_bytes,
            open_marks,
            meter_report_marks,
            idempotency_marks,
            control_marks,
        })
    }
}

#[cfg(test)]
fn push_bounded_event(
    buffered: &mut Vec<BufferedLedgerEvent>,
    encoded_bytes: &mut usize,
    event: BufferedLedgerEvent,
) -> Result<bool, SessionError> {
    let event_bytes = event.encoded_len()?;
    if event_bytes > MAX_FLUSH_ENCODED_BYTES {
        return Err(SessionError::LimitExceeded {
            resource: "flush_row_encoded_bytes",
            limit: MAX_FLUSH_ENCODED_BYTES,
            observed_at_least: event_bytes,
        });
    }
    let next_bytes = encoded_bytes
        .checked_add(event_bytes)
        .ok_or(SessionError::LimitExceeded {
            resource: "flush_encoded_bytes",
            limit: MAX_FLUSH_ENCODED_BYTES,
            observed_at_least: usize::MAX,
        })?;
    if buffered.len() == MAX_FLUSH_ROWS || next_bytes > MAX_FLUSH_ENCODED_BYTES {
        return Ok(false);
    }
    buffered.push(event);
    *encoded_bytes = next_bytes;
    Ok(true)
}

#[derive(Default)]
struct ScopeState {
    sessions: BTreeSet<u64>,
    dirty_open_receipts: BTreeSet<SessionOpenId>,
    dirty_causal: BTreeSet<(u64, DirtyCausalMutation)>,
    dirty_submission_failures: BTreeSet<(u64, SubmissionRequestId)>,
    dirty_control: BTreeSet<(i64, DirtyControlMutation)>,
    // Immutable evidence-bearing rows make page/flush snapshots shallow while
    // the governor lock is held. Public cloning and JSON materialization occur
    // only after that lock is released.
    events: Vec<Arc<DegradationEvent>>,
    flushed_events: usize,
    sink: Option<fs_ledger::LedgerInstanceId>,
    flush_generation: i64,
    in_flight: Option<u64>,
    revision: u64,
    next_flush_lane: u8,
    reserved_pause_completions: usize,
    retained_bytes: usize,
}

fn has_dirty_submission_predecessor(
    scope: &ScopeState,
    pause_request: PauseRequestId,
    selected_submissions: &BTreeSet<SubmissionRequestId>,
) -> bool {
    let is_draining_generation = |request_id: SubmissionRequestId| {
        request_id.session == pause_request.session
            && request_id.generation == pause_request.gate_generation
    };
    let is_selected = |request_id: SubmissionRequestId| selected_submissions.contains(&request_id);
    scope.dirty_causal.iter().any(|(_, mutation)| {
        matches!(
            mutation,
            DirtyCausalMutation::Submission(request_id)
                if is_draining_generation(*request_id) && !is_selected(*request_id)
        )
    }) || scope
        .dirty_submission_failures
        .iter()
        .any(|(_, request_id)| is_draining_generation(*request_id) && !is_selected(*request_id))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PendingPause {
    request_id: PauseRequestId,
    pressure_action_id: PressureActionId,
    reserved_retained_bytes: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CompletedPause {
    request_id: PauseRequestId,
    checkpoint_byte_len: usize,
    checkpoint_digest: fs_blake3::ContentHash,
    completion_event_index: usize,
    completion_ordinal: i64,
    resume_generation: u64,
    gate_binding: fs_blake3::ContentHash,
    acknowledgement_hash: fs_blake3::ContentHash,
}

#[derive(Debug, Clone, Copy)]
struct PauseAcknowledgementReplay {
    completion_event_index: usize,
    resume_generation: u64,
    gate_binding: fs_blake3::ContentHash,
    content_hash: fs_blake3::ContentHash,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GatePhase {
    Running,
    ReadyToResume,
}

#[derive(Default)]
struct Inner {
    open_requests: BTreeMap<SessionOpenId, OpenReplay>,
    session_open_ids: BTreeMap<u64, SessionOpenId>,
    tokens: BTreeMap<u64, Arc<CapabilityToken>>,
    /// Session-OWNED cancellation gates, bound at open (gp3.13): the
    /// only route to a pause request, so a foreign gate is
    /// unrepresentable at the pressure API.
    gates: BTreeMap<u64, Arc<CancelGate>>,
    /// Current cancellation-gate generation for every gated session.
    gate_generations: BTreeMap<u64, u64>,
    /// Whether the current gate is running or awaiting explicit activation.
    gate_phases: BTreeMap<u64, GatePhase>,
    /// Pause requests awaiting a checkpoint acknowledgement, keyed by
    /// session → ordinal of the Requested event.
    pending_pause: BTreeMap<u64, PendingPause>,
    /// Last completed request per session for idempotent response replay.
    completed_pause: BTreeMap<u64, CompletedPause>,
    pause_acknowledgements: BTreeMap<PauseRequestId, PauseAcknowledgementReplay>,
    resume_activations: BTreeMap<ResumeActivationId, ResumeActivationReceipt>,
    meters: BTreeMap<u64, SessionMeters>,
    meter_reports: BTreeMap<MeterReportId, MeterReceipt>,
    meter_report_ids: BTreeMap<u64, BTreeSet<MeterReportId>>,
    reserved_meter_reports: BTreeMap<u64, usize>,
    pressure_actions: BTreeMap<PressureActionId, PressureReplay>,
    pressure_action_ids: BTreeMap<u64, BTreeSet<PressureActionId>>,
    idempotency: BTreeMap<SubmissionRequestId, IdemState>,
    idempotency_keys: BTreeMap<u64, BTreeMap<fs_blake3::ContentHash, SubmissionRequestId>>,
    submission_admissions: BTreeMap<u64, SubmissionRequestId>,
    pending_submissions: BTreeMap<u64, usize>,
    scopes: BTreeMap<String, ScopeState>,
    next_submission_ordinal: u64,
    next_meter_commit_ordinal: u64,
    reserved_meter_ordinals: usize,
    next_ordinal: i64,
    reserved_pause_ordinals: usize,
    next_flush_reservation: u64,
    retained_bytes: usize,
    durable_claims_at_open: u64,
    recovered_durable_claims: BTreeSet<fs_blake3::ContentHash>,
    program_risk_publications: BTreeSet<fs_blake3::ContentHash>,
}

fn checked_retained_add(current: usize, added: usize) -> usize {
    current.saturating_add(added)
}

fn ensure_retained_capacity(
    inner: &Inner,
    ledger_scope: &str,
    added: usize,
) -> Result<(), SessionError> {
    let scope_current = inner
        .scopes
        .get(ledger_scope)
        .map_or(0, |scope| scope.retained_bytes);
    let scope_next = checked_retained_add(scope_current, added);
    if scope_next > MAX_RETAINED_BYTES_PER_SCOPE {
        return Err(SessionError::LimitExceeded {
            resource: "retained_bytes_per_scope",
            limit: MAX_RETAINED_BYTES_PER_SCOPE,
            observed_at_least: scope_next,
        });
    }
    let governor_next = checked_retained_add(inner.retained_bytes, added);
    if governor_next > MAX_RETAINED_BYTES_PER_GOVERNOR {
        return Err(SessionError::LimitExceeded {
            resource: "retained_bytes_per_governor",
            limit: MAX_RETAINED_BYTES_PER_GOVERNOR,
            observed_at_least: governor_next,
        });
    }
    Ok(())
}

fn commit_retained_bytes(inner: &mut Inner, ledger_scope: &str, added: usize) {
    inner.retained_bytes = inner
        .retained_bytes
        .checked_add(added)
        .expect("retained-capacity preflight prevents governor overflow");
    let scope = inner
        .scopes
        .get_mut(ledger_scope)
        .expect("registered session scope");
    scope.retained_bytes = scope
        .retained_bytes
        .checked_add(added)
        .expect("retained-capacity preflight prevents scope overflow");
}

fn release_retained_bytes(inner: &mut Inner, ledger_scope: &str, released: usize) {
    inner.retained_bytes = inner
        .retained_bytes
        .checked_sub(released)
        .expect("released bytes were previously reserved");
    let scope = inner
        .scopes
        .get_mut(ledger_scope)
        .expect("registered session scope");
    scope.retained_bytes = scope
        .retained_bytes
        .checked_sub(released)
        .expect("released scope bytes were previously reserved");
}

fn bump_scope_revision(inner: &mut Inner, ledger_scope: &str) {
    let scope = inner
        .scopes
        .get_mut(ledger_scope)
        .expect("registered session scope");
    // A saturated revision makes `remaining_dirty` conservatively stay true;
    // collection and ordinal bounds are reached vastly earlier.
    scope.revision = scope.revision.saturating_add(1);
}

/// The governor. `Send + Sync`: hot paths are mutex-guarded in-memory
/// state; ledger persistence is the explicit single-threaded flush.
pub struct Governor {
    id: fs_blake3::ContentHash,
    durable_sink: Option<fs_ledger::LedgerInstanceId>,
    inner: Mutex<Inner>,
    #[cfg(test)]
    materialization_barrier: Mutex<Option<Arc<std::sync::Barrier>>>,
}

/// Process-local exclusion for one program-risk singleton publication.
///
/// The ledger terminal remains the durable arbitration point. This guard
/// prevents two callers sharing the same live governor from materializing
/// duplicate lineage operations before either reaches that terminal.
pub(crate) struct ProgramRiskPublicationGuard<'a> {
    governor: &'a Governor,
    authority: fs_blake3::ContentHash,
}

impl Drop for ProgramRiskPublicationGuard<'_> {
    fn drop(&mut self) {
        self.governor
            .inner
            .lock()
            .expect("governor lock")
            .program_risk_publications
            .remove(&self.authority);
    }
}

impl Default for Governor {
    fn default() -> Self {
        Governor::new()
    }
}

impl Governor {
    /// An empty governor.
    #[must_use]
    pub fn new() -> Self {
        Governor {
            id: ephemeral_governor_id(),
            durable_sink: None,
            inner: Mutex::new(Inner::default()),
            #[cfg(test)]
            materialization_barrier: Mutex::new(None),
        }
    }

    /// Construct a restart-stable governor identity bound to one physical
    /// ledger instance and one explicit caller-persisted nonce.
    ///
    /// Repeating this call after reopening the same ledger with the same nonce
    /// reconstructs the exact authority namespace. The constructor snapshots
    /// its complete durable claim count; fresh mutation stays fenced until the
    /// typed recovery APIs have rebuilt that governor-wide history. A
    /// replacement or foreign ledger derives a different identity.
    ///
    /// # Errors
    /// Corrupt or unavailable physical-ledger identity fails closed.
    pub fn new_durable(
        ledger: &fs_ledger::Ledger,
        nonce: DurableGovernorNonce,
    ) -> Result<Self, SessionError> {
        let sink = ledger
            .checked_instance_id()
            .map_err(|error| SessionError::Persistence {
                what: format!("durable governor ledger identity validation failed: {error}"),
            })?;
        let mut payload = Vec::new();
        payload.extend_from_slice(&sink.as_bytes());
        payload.extend_from_slice(&nonce.as_bytes());
        let id = fs_blake3::hash_domain(DURABLE_GOVERNOR_ID_DOMAIN, &payload);
        let durable_claims_at_open =
            ledger
                .session_mutation_claim_count(id)
                .map_err(|error| SessionError::Persistence {
                    what: format!("durable governor recovery-fence snapshot failed: {error}"),
                })?;
        let inner = Inner {
            durable_claims_at_open,
            ..Inner::default()
        };
        Ok(Self {
            id,
            durable_sink: Some(sink),
            inner: Mutex::new(inner),
            #[cfg(test)]
            materialization_barrier: Mutex::new(None),
        })
    }

    #[cfg(test)]
    fn set_materialization_barrier(&self, barrier: Option<Arc<std::sync::Barrier>>) {
        *self
            .materialization_barrier
            .lock()
            .expect("materialization barrier lock") = barrier;
    }

    #[cfg(test)]
    fn wait_at_materialization_barrier_for_test(&self) {
        let barrier = self
            .materialization_barrier
            .lock()
            .expect("materialization barrier lock")
            .clone();
        if let Some(barrier) = barrier {
            barrier.wait();
        }
    }

    #[cfg(not(test))]
    fn wait_at_materialization_barrier_for_test(&self) {}

    fn ensure_durable_recovery_complete(&self, inner: &Inner) -> Result<(), SessionError> {
        if self.durable_sink.is_none() {
            return Ok(());
        }
        let recovered = u64::try_from(inner.recovered_durable_claims.len()).map_err(|_| {
            SessionError::Persistence {
                what: "recovered durable-claim count overflowed u64".to_string(),
            }
        })?;
        if recovered < inner.durable_claims_at_open {
            return Err(SessionError::DurableRecoveryIncomplete {
                remaining_claims: inner.durable_claims_at_open - recovered,
            });
        }
        Ok(())
    }

    fn mark_durable_claim_recovered(&self, inner: &mut Inner, authority: fs_blake3::ContentHash) {
        if self.durable_sink.is_some() {
            inner.recovered_durable_claims.insert(authority);
        }
    }

    /// Opaque governor identity carried by every typed mutation authority.
    #[must_use]
    pub const fn identity(&self) -> fs_blake3::ContentHash {
        self.id
    }

    /// Mint bounded retry authority before opening a session.
    ///
    /// # Errors
    /// Blank or oversized client keys are refused without retaining state.
    pub fn session_open_id(
        &self,
        session: SessionId,
        client_key: &str,
    ) -> Result<SessionOpenId, SessionError> {
        let key_hash =
            bounded_request_digest("session_open_key_bytes", SESSION_OPEN_ID_DOMAIN, client_key)?;
        let mut payload = Vec::new();
        push_hash(&mut payload, self.id);
        payload.extend_from_slice(&session.0.to_le_bytes());
        push_hash(&mut payload, key_hash);
        Ok(SessionOpenId {
            governor_id: self.id,
            session,
            content_hash: fs_blake3::hash_domain(SESSION_OPEN_ID_DOMAIN, &payload),
        })
    }

    #[allow(clippy::too_many_lines)] // One rollback-free authority registration transaction.
    fn register_session(
        &self,
        open_id: SessionOpenId,
        mut token: CapabilityToken,
        gate: Option<Arc<CancelGate>>,
        recovering: bool,
    ) -> Result<SessionOpenReceipt, SessionError> {
        if open_id.governor_id != self.id || open_id.session != token.session {
            return Err(SessionError::MutationAuthorityMismatch {
                kind: "session-open",
                id: open_id.content_hash,
            });
        }
        token.validate_operator_grants()?;
        CapabilityToken::validate_ledger_scope(&token.ledger_scope)?;
        validate_resource("core-seconds grant", token.core_s)?;
        validate_resource("wall-seconds grant", token.wall_s)?;
        // Public tokens are caller-constructed, so valid short strings and
        // vectors may still carry attacker-chosen spare capacities. Rebuild
        // them from validated slices before retention; allocator rounding is
        // then bounded by the admitted content instead of caller history.
        token.ops = token
            .ops
            .iter()
            .map(|grant| grant.as_str().to_owned())
            .collect();
        let caller_scope = std::mem::take(&mut token.ledger_scope);
        token.ledger_scope = String::from(caller_scope.as_str());
        let session = token.session.0;
        let ledger_scope = token.ledger_scope.clone();
        let token_digest = capability_token_identity(&token);
        let mut g = self.inner.lock().expect("governor lock");
        if let Some(replay) = g.open_requests.get(&open_id) {
            let same_gate = match (&replay.gate, &gate) {
                (None, None) => true,
                (Some(left), Some(right)) => Arc::ptr_eq(left, right),
                (None, Some(_)) | (Some(_), None) => false,
            };
            if replay.token_digest == token_digest && same_gate {
                return Ok(replay.receipt.as_ref().clone());
            }
            return Err(SessionError::MutationConflict {
                kind: "session-open",
                id: open_id.content_hash,
            });
        }
        if !recovering {
            self.ensure_durable_recovery_complete(&g)?;
        }
        if gate.as_ref().is_some_and(|gate| gate.is_requested()) {
            return Err(SessionError::PreRequestedGate { id: session });
        }
        if g.tokens.contains_key(&session) {
            return Err(SessionError::SessionAlreadyOpen { id: session });
        }
        if g.tokens.len() >= MAX_SESSIONS_PER_GOVERNOR {
            return Err(SessionError::LimitExceeded {
                resource: "sessions_per_governor",
                limit: MAX_SESSIONS_PER_GOVERNOR,
                observed_at_least: g.tokens.len().saturating_add(1),
            });
        }
        let scope_session_count = g
            .scopes
            .get(&ledger_scope)
            .map_or(0, |scope| scope.sessions.len());
        if scope_session_count >= MAX_SESSIONS_PER_SCOPE {
            return Err(SessionError::LimitExceeded {
                resource: "sessions_per_scope",
                limit: MAX_SESSIONS_PER_SCOPE,
                observed_at_least: scope_session_count.saturating_add(1),
            });
        }
        let operator_bytes: usize = token.ops.iter().map(String::len).sum();
        let retained_bytes = ledger_scope
            .len()
            .checked_add(operator_bytes)
            .and_then(|bytes| {
                bytes.checked_add(if g.scopes.contains_key(&ledger_scope) {
                    0
                } else {
                    ledger_scope.len()
                })
            })
            .and_then(|bytes| bytes.checked_add(OPEN_REQUEST_RETAINED_BYTES))
            .ok_or(SessionError::LimitExceeded {
                resource: "retained_bytes_per_governor",
                limit: MAX_RETAINED_BYTES_PER_GOVERNOR,
                observed_at_least: usize::MAX,
            })?;
        ensure_retained_capacity(&g, &ledger_scope, retained_bytes)?;
        let next_revision = g
            .scopes
            .get(&ledger_scope)
            .map_or(1, |scope| scope.revision.saturating_add(1));
        let gate_identity = gate.as_ref().map(|_| session_gate_binding(open_id));
        let permit = ScopeFlushPermit {
            governor_id: self.id,
            ledger_scope: ledger_scope.clone(),
        };
        let receipt = SessionOpenReceipt {
            open_id,
            token_digest,
            gate_identity,
            permit,
            content_hash: session_open_receipt_hash(
                open_id,
                token_digest,
                gate_identity,
                &ledger_scope,
            ),
        };
        let token = Arc::new(token);
        let retained_receipt = Arc::new(receipt.clone());
        g.meters.insert(session, SessionMeters::default());
        g.meter_report_ids.insert(session, BTreeSet::new());
        g.reserved_meter_reports.insert(session, 0);
        g.pressure_action_ids.insert(session, BTreeSet::new());
        g.idempotency_keys.insert(session, BTreeMap::new());
        g.pending_submissions.insert(session, 0);
        g.session_open_ids.insert(session, open_id);
        g.tokens.insert(session, token);
        if let Some(bound_gate) = &gate {
            g.gates.insert(session, Arc::clone(bound_gate));
            g.gate_generations.insert(session, 0);
            g.gate_phases.insert(session, GatePhase::Running);
        }
        let scope = g.scopes.entry(ledger_scope.clone()).or_default();
        scope.sessions.insert(session);
        scope.dirty_open_receipts.insert(open_id);
        scope.revision = next_revision;
        g.open_requests.insert(
            open_id,
            OpenReplay {
                token_digest,
                gate,
                receipt: retained_receipt,
            },
        );
        commit_retained_bytes(&mut g, &ledger_scope, retained_bytes);
        Ok(receipt)
    }

    /// Register a session's token (issuance). Session ids are single-use for
    /// the lifetime of this governor; duplicate registration fails closed.
    ///
    /// # Errors
    /// - [`SessionError::InvalidLedgerScope`] when the token's namespace is not
    ///   canonical and bounded.
    /// - [`SessionError::InvalidResource`] when a floating-point time grant is
    ///   not finite and non-negative.
    /// - [`SessionError::SessionAlreadyOpen`] when the id is already
    ///   registered.
    /// - [`SessionError::LimitExceeded`] when the governor-wide or scoped
    ///   session cap has been reached.
    ///
    /// Integer memory/core grants are structurally bounded. Rejection happens
    /// before any session state is mutated.
    pub fn open_session(
        &self,
        open_id: SessionOpenId,
        token: CapabilityToken,
    ) -> Result<SessionOpenReceipt, SessionError> {
        self.register_session(open_id, token, None, false)
    }

    /// Register a session's token WITH its cancellation capability
    /// (bead gp3.13): the gate is owned by the governor from open, and
    /// level-3 memory pressure resolves it by `SessionId` — passing
    /// someone else's gate to a pressure action is unrepresentable.
    /// Sessions opened without a gate refuse level-3 pressure.
    ///
    /// # Errors
    /// The same [`SessionError::InvalidLedgerScope`],
    /// [`SessionError::InvalidResource`],
    /// [`SessionError::SessionAlreadyOpen`], and
    /// [`SessionError::LimitExceeded`] refusals as
    /// [`Governor::open_session`].
    pub fn open_session_gated(
        &self,
        open_id: SessionOpenId,
        token: CapabilityToken,
        gate: Arc<CancelGate>,
    ) -> Result<SessionOpenReceipt, SessionError> {
        self.register_session(open_id, token, Some(gate), false)
    }

    /// The token for a session.
    ///
    /// # Errors
    /// [`SessionError::UnknownSession`].
    pub fn token(&self, session: SessionId) -> Result<CapabilityToken, SessionError> {
        self.inner
            .lock()
            .expect("governor lock")
            .tokens
            .get(&session.0)
            .map(|token| token.as_ref().clone())
            .ok_or(SessionError::UnknownSession { id: session.0 })
    }

    pub(crate) fn ensure_program_risk_mutation_ready(&self) -> Result<(), SessionError> {
        let inner = self.inner.lock().expect("governor lock");
        self.ensure_durable_recovery_complete(&inner)
    }

    pub(crate) fn reserve_program_risk_publication(
        &self,
        authority: fs_blake3::ContentHash,
    ) -> Result<ProgramRiskPublicationGuard<'_>, SessionError> {
        let mut inner = self.inner.lock().expect("governor lock");
        if !inner.program_risk_publications.insert(authority) {
            return Err(SessionError::ProgramRiskReportInFlight { id: authority });
        }
        Ok(ProgramRiskPublicationGuard {
            governor: self,
            authority,
        })
    }

    pub(crate) fn program_risk_session_context(
        &self,
        receipt: &SessionOpenReceipt,
        attempted_sink: fs_ledger::LedgerInstanceId,
    ) -> Result<(SessionId, String, u64), SessionError> {
        let open_id = receipt.open_id;
        if open_id.governor_id != self.id {
            return Err(SessionError::MutationAuthorityMismatch {
                kind: recovery::KIND_PROGRAM_RISK_REPORT,
                id: open_id.content_hash,
            });
        }
        let inner = self.inner.lock().expect("governor lock");
        let known_open =
            inner
                .session_open_ids
                .get(&open_id.session.0)
                .ok_or(SessionError::UnknownSession {
                    id: open_id.session.0,
                })?;
        let replay =
            inner
                .open_requests
                .get(known_open)
                .ok_or_else(|| SessionError::Persistence {
                    what: format!(
                        "session {} lost its immutable open receipt",
                        open_id.session.0
                    ),
                })?;
        if *known_open != open_id || replay.receipt.as_ref() != receipt {
            return Err(SessionError::MutationAuthorityMismatch {
                kind: recovery::KIND_PROGRAM_RISK_REPORT,
                id: open_id.content_hash,
            });
        }
        let pending_submissions = inner
            .pending_submissions
            .get(&open_id.session.0)
            .copied()
            .unwrap_or(0);
        let pause_pending = inner.pending_pause.contains_key(&open_id.session.0);
        if pending_submissions != 0 || pause_pending {
            return Err(SessionError::SessionNotQuiescent {
                id: open_id.session.0,
                pending_submissions,
                pause_pending,
            });
        }
        let token = inner
            .tokens
            .get(&open_id.session.0)
            .expect("validated open receipt retains its token");
        let scope = inner
            .scopes
            .get(&token.ledger_scope)
            .expect("validated token retains its scope");
        match scope.sink {
            Some(bound_sink) if bound_sink != attempted_sink => {
                return Err(SessionError::LedgerScopeSinkMismatch {
                    scope: token.ledger_scope.clone(),
                    bound_sink,
                    attempted_sink,
                });
            }
            None => {
                return Err(SessionError::Persistence {
                    what: format!(
                        "session {} scope {:?} must be flushed to the target ledger before its program-risk session-end report",
                        open_id.session.0, token.ledger_scope
                    ),
                });
            }
            Some(_) => {}
        }
        let generation = inner
            .gate_generations
            .get(&open_id.session.0)
            .copied()
            .unwrap_or(0);
        Ok((open_id.session, token.ledger_scope.clone(), generation))
    }

    pub(crate) fn mark_program_risk_report_recovered(&self, authority: fs_blake3::ContentHash) {
        let mut inner = self.inner.lock().expect("governor lock");
        self.mark_durable_claim_recovered(&mut inner, authority);
    }

    fn mutation_context(
        &self,
        session: SessionId,
    ) -> Result<(fs_blake3::ContentHash, u64), SessionError> {
        let g = self.inner.lock().expect("governor lock");
        let session_open = Self::current_open_identity(&g, session)?;
        let generation = g.gate_generations.get(&session.0).copied().unwrap_or(0);
        Ok((session_open, generation))
    }

    fn current_open_identity(
        g: &Inner,
        session: SessionId,
    ) -> Result<fs_blake3::ContentHash, SessionError> {
        let open_id = g
            .session_open_ids
            .get(&session.0)
            .ok_or(SessionError::UnknownSession { id: session.0 })?;
        g.open_requests
            .get(open_id)
            .map(|replay| replay.receipt.content_hash)
            .ok_or_else(|| SessionError::Persistence {
                what: format!(
                    "session {} lost its immutable open receipt identity",
                    session.0
                ),
            })
    }

    /// Mint a bounded authority for one exact-bit meter report.
    ///
    /// # Errors
    /// Blank/oversized keys, unknown sessions, or corrupt immutable open state
    /// are refused without minting an authority.
    pub fn meter_report_id(
        &self,
        session: SessionId,
        client_key: &str,
    ) -> Result<MeterReportId, SessionError> {
        let key_hash =
            bounded_request_digest("meter_report_key_bytes", METER_REPORT_ID_DOMAIN, client_key)?;
        let (session_open, generation) = self.mutation_context(session)?;
        let mut payload = Vec::new();
        push_hash(&mut payload, self.id);
        payload.extend_from_slice(&session.0.to_le_bytes());
        push_hash(&mut payload, session_open);
        payload.extend_from_slice(&generation.to_le_bytes());
        push_hash(&mut payload, key_hash);
        Ok(MeterReportId {
            governor_id: self.id,
            session,
            session_open,
            generation,
            content_hash: fs_blake3::hash_domain(METER_REPORT_ID_DOMAIN, &payload),
        })
    }

    fn validate_meter_authority(
        &self,
        g: &Inner,
        report_id: MeterReportId,
    ) -> Result<Arc<CapabilityToken>, SessionError> {
        let token =
            g.tokens
                .get(&report_id.session.0)
                .cloned()
                .ok_or(SessionError::UnknownSession {
                    id: report_id.session.0,
                })?;
        let current_open = Self::current_open_identity(g, report_id.session)?;
        if report_id.governor_id != self.id || report_id.session_open != current_open {
            return Err(SessionError::MutationAuthorityMismatch {
                kind: "meter-report",
                id: report_id.content_hash,
            });
        }
        let current_generation = g
            .gate_generations
            .get(&report_id.session.0)
            .copied()
            .unwrap_or(0);
        if report_id.generation != current_generation {
            return Err(SessionError::StaleMutationGeneration {
                kind: "meter-report",
                id: report_id.session.0,
                supplied: report_id.generation,
                current: current_generation,
            });
        }
        Ok(token)
    }

    #[allow(clippy::too_many_lines)] // Validation, reservation, transition, and receipt commit are one atomic path.
    fn commit_meter_locked(
        &self,
        g: &mut Inner,
        report_id: MeterReportId,
        delta: Charge,
        consumes_reservation: bool,
    ) -> Result<MeterReceipt, SessionError> {
        if report_id.governor_id != self.id {
            return Err(SessionError::MutationAuthorityMismatch {
                kind: "meter-report",
                id: report_id.content_hash,
            });
        }
        if let Some(receipt) = g.meter_reports.get(&report_id) {
            if same_charge(receipt.delta, delta) {
                return Ok(receipt.clone());
            }
            return Err(SessionError::MutationConflict {
                kind: "meter-report",
                id: report_id.content_hash,
            });
        }
        self.ensure_durable_recovery_complete(g)?;
        let token = self.validate_meter_authority(g, report_id)?;
        let report_count = g
            .meter_report_ids
            .get(&report_id.session.0)
            .map_or(0, BTreeSet::len);
        let reserved = g
            .reserved_meter_reports
            .get(&report_id.session.0)
            .copied()
            .unwrap_or(0);
        let occupied = report_count
            .checked_add(reserved)
            .ok_or(SessionError::LimitExceeded {
                resource: "meter_reports_per_session",
                limit: MAX_METER_REPORTS_PER_SESSION,
                observed_at_least: usize::MAX,
            })?;
        if (!consumes_reservation && occupied >= MAX_METER_REPORTS_PER_SESSION)
            || (consumes_reservation && reserved == 0)
        {
            return Err(SessionError::LimitExceeded {
                resource: "meter_reports_per_session",
                limit: MAX_METER_REPORTS_PER_SESSION,
                observed_at_least: occupied.saturating_add(usize::from(!consumes_reservation)),
            });
        }
        if consumes_reservation && g.reserved_meter_ordinals == 0 {
            return Err(SessionError::Persistence {
                what: "a pending submission lost its reserved meter-commit ordinal capacity"
                    .to_string(),
            });
        }
        let unowned_advance = usize::from(!consumes_reservation);
        let reserved_advance = g
            .reserved_meter_ordinals
            .checked_add(unowned_advance)
            .and_then(|advance| u64::try_from(advance).ok())
            .ok_or(SessionError::LimitExceeded {
                resource: "meter_commit_ordinal",
                limit: i64::MAX as usize,
                observed_at_least: usize::MAX,
            })?;
        if g.next_meter_commit_ordinal
            .checked_add(reserved_advance)
            .is_none_or(|last| last > i64::MAX as u64)
        {
            return Err(SessionError::LimitExceeded {
                resource: "meter_commit_ordinal",
                limit: i64::MAX as usize,
                observed_at_least: i64::MAX as usize,
            });
        }
        if !consumes_reservation {
            ensure_retained_capacity(g, &token.ledger_scope, MAX_METER_RECEIPT_RETAINED_BYTES)?;
        }
        let before = g
            .meters
            .get(&report_id.session.0)
            .cloned()
            .unwrap_or_default();
        let (next, enforcement) = meter_transition(&token, &before, delta)?;
        let commit_ordinal =
            g.next_meter_commit_ordinal
                .checked_add(1)
                .ok_or(SessionError::LimitExceeded {
                    resource: "meter_commit_ordinal",
                    limit: i64::MAX as usize,
                    observed_at_least: usize::MAX,
                })?;
        if commit_ordinal > i64::MAX as u64 {
            return Err(SessionError::LimitExceeded {
                resource: "meter_commit_ordinal",
                limit: i64::MAX as usize,
                observed_at_least: i64::MAX as usize,
            });
        }
        let before_snapshot = before.snapshot();
        let after_snapshot = next.snapshot();
        let receipt = MeterReceipt {
            report_id,
            commit_ordinal,
            delta,
            before: before_snapshot,
            after: after_snapshot,
            content_hash: meter_receipt_hash(
                report_id,
                commit_ordinal,
                delta,
                before_snapshot,
                after_snapshot,
                &enforcement,
            ),
            enforcement,
        };
        g.next_meter_commit_ordinal = commit_ordinal;
        if consumes_reservation {
            *g.reserved_meter_reports
                .get_mut(&report_id.session.0)
                .expect("open session has meter reservation count") -= 1;
            g.reserved_meter_ordinals -= 1;
        } else {
            commit_retained_bytes(g, &token.ledger_scope, MAX_METER_RECEIPT_RETAINED_BYTES);
        }
        g.meters.insert(report_id.session.0, next);
        g.meter_report_ids
            .get_mut(&report_id.session.0)
            .expect("open session has meter-report index")
            .insert(report_id);
        g.meter_reports.insert(report_id, receipt.clone());
        let scope = g
            .scopes
            .get_mut(&token.ledger_scope)
            .expect("registered session scope");
        if !consumes_reservation {
            scope
                .dirty_causal
                .insert((commit_ordinal, DirtyCausalMutation::Meter(report_id)));
        }
        bump_scope_revision(g, &token.ledger_scope);
        Ok(receipt)
    }

    /// Commit or exactly replay one metering report. The payload comparison is
    /// exact-bit and a duplicate changes no meter, counter, ordinal, or cursor.
    ///
    /// # Errors
    /// Foreign/stale/conflicting authority, invalid charge, capacity, and
    /// corrupt session-state failures are returned before partial mutation.
    pub fn charge(
        &self,
        report_id: MeterReportId,
        delta: Charge,
    ) -> Result<MeterReceipt, SessionError> {
        let mut g = self.inner.lock().expect("governor lock");
        self.commit_meter_locked(&mut g, report_id, delta, false)
    }

    /// Mint a retry authority from a stable caller key and canonical program.
    /// The caller key selects one mutation slot; changing the program under
    /// that slot is detected as a conflict when submitted.
    ///
    /// # Errors
    /// Blank/oversized inputs, unknown sessions, or corrupt immutable open
    /// state are refused without minting an authority.
    pub fn submission_request_id(
        &self,
        session: SessionId,
        agent_key: &str,
        program_text: &str,
    ) -> Result<SubmissionRequestId, SessionError> {
        let key_hash = bounded_request_digest(
            "idempotency_agent_key_bytes",
            IDEMPOTENCY_AGENT_DOMAIN,
            agent_key,
        )?;
        let request_hash = bounded_request_digest(
            "idempotency_program_text_bytes",
            IDEMPOTENCY_PROGRAM_DOMAIN,
            program_text,
        )?;
        let (session_open, generation) = self.mutation_context(session)?;
        let mut payload = Vec::new();
        push_hash(&mut payload, self.id);
        payload.extend_from_slice(&session.0.to_le_bytes());
        push_hash(&mut payload, session_open);
        payload.extend_from_slice(&generation.to_le_bytes());
        push_hash(&mut payload, key_hash);
        Ok(SubmissionRequestId {
            governor_id: self.id,
            session,
            session_open,
            generation,
            key_hash,
            request_hash,
            content_hash: fs_blake3::hash_domain(SUBMISSION_REQUEST_ID_DOMAIN, &payload),
        })
    }

    fn submission_meter_report_id(request_id: SubmissionRequestId) -> MeterReportId {
        let mut payload = Vec::new();
        push_hash(&mut payload, request_id.content_hash);
        MeterReportId {
            governor_id: request_id.governor_id,
            session: request_id.session,
            session_open: request_id.session_open,
            generation: request_id.generation,
            content_hash: fs_blake3::hash_domain(METER_REPORT_ID_DOMAIN, &payload),
        }
    }

    #[allow(clippy::too_many_lines)] // Preflight every accounting invariant before one exact rollback.
    fn rollback_submission_admission(
        &self,
        request_id: SubmissionRequestId,
        admission_ordinal: u64,
        ledger_scope: &str,
    ) -> Result<(), SessionError> {
        let session = request_id.session;
        let mut g = self.inner.lock().expect("governor lock");
        let (reserved_terminal_bytes, reserved_meter_bytes) = match g.idempotency.get(&request_id) {
            Some(IdemState::Pending {
                admission_ordinal: pending_ordinal,
                request_id: pending_request,
                reserved_terminal_bytes,
                reserved_meter_bytes,
                durable_permit: None,
            }) if *pending_ordinal == admission_ordinal && *pending_request == request_id => {
                (*reserved_terminal_bytes, *reserved_meter_bytes)
            }
            _ => {
                return Err(SessionError::Persistence {
                    what: format!(
                        "session {} cannot roll back submission admission {} because its exact unclaimed Pending state changed",
                        session.0, admission_ordinal
                    ),
                });
            }
        };
        let retained_bytes = SUBMISSION_REQUEST_RETAINED_BYTES
            .checked_add(reserved_terminal_bytes)
            .and_then(|bytes| bytes.checked_add(reserved_meter_bytes))
            .ok_or(SessionError::LimitExceeded {
                resource: "retained_bytes_per_scope",
                limit: MAX_RETAINED_BYTES_PER_SCOPE,
                observed_at_least: usize::MAX,
            })?;
        let pending_count = g
            .pending_submissions
            .get(&session.0)
            .copied()
            .ok_or_else(|| SessionError::Persistence {
                what: format!(
                    "session {} lost its pending-submission counter during claim rollback",
                    session.0
                ),
            })?;
        let reserved_reports = g
            .reserved_meter_reports
            .get(&session.0)
            .copied()
            .ok_or_else(|| SessionError::Persistence {
                what: format!(
                    "session {} lost its meter-reservation counter during claim rollback",
                    session.0
                ),
            })?;
        if pending_count == 0 || reserved_reports == 0 || g.reserved_meter_ordinals == 0 {
            return Err(SessionError::Persistence {
                what: format!(
                    "session {} has exhausted reservation counters during claim rollback",
                    session.0
                ),
            });
        }
        let key_index =
            g.idempotency_keys
                .get(&session.0)
                .ok_or_else(|| SessionError::Persistence {
                    what: format!(
                        "session {} lost its submission-key index during claim rollback",
                        session.0
                    ),
                })?;
        if key_index.get(&request_id.key_hash) != Some(&request_id) {
            return Err(SessionError::Persistence {
                what: format!(
                    "session {} submission key changed during claim rollback",
                    session.0
                ),
            });
        }
        if g.submission_admissions.get(&admission_ordinal) != Some(&request_id) {
            return Err(SessionError::Persistence {
                what: format!(
                    "submission admission ordinal {admission_ordinal} changed owner during claim rollback"
                ),
            });
        }
        let scope_retained = g
            .scopes
            .get(ledger_scope)
            .map(|scope| scope.retained_bytes)
            .ok_or_else(|| SessionError::Persistence {
                what: format!("scope {ledger_scope} disappeared during claim rollback"),
            })?;
        if scope_retained < retained_bytes || g.retained_bytes < retained_bytes {
            return Err(SessionError::Persistence {
                what: format!(
                    "scope {ledger_scope} retained-byte accounting underflow during claim rollback"
                ),
            });
        }

        g.idempotency.remove(&request_id);
        g.idempotency_keys
            .get_mut(&session.0)
            .expect("submission-key index checked above")
            .remove(&request_id.key_hash);
        g.submission_admissions.remove(&admission_ordinal);
        *g.pending_submissions
            .get_mut(&session.0)
            .expect("pending-submission counter checked above") -= 1;
        *g.reserved_meter_reports
            .get_mut(&session.0)
            .expect("meter-reservation counter checked above") -= 1;
        g.reserved_meter_ordinals -= 1;
        if g.next_submission_ordinal == admission_ordinal {
            g.next_submission_ordinal -= 1;
        }
        release_retained_bytes(&mut g, ledger_scope, retained_bytes);
        Ok(())
    }

    fn attach_submission_permit(
        &self,
        request_id: SubmissionRequestId,
        admission_ordinal: u64,
        permit: fs_ledger::session_registry::SessionClaimPermit,
    ) -> Result<(), SessionError> {
        let mut g = self.inner.lock().expect("governor lock");
        match g.idempotency.get_mut(&request_id) {
            Some(IdemState::Pending {
                admission_ordinal: pending_ordinal,
                request_id: pending_request,
                durable_permit,
                ..
            }) if *pending_ordinal == admission_ordinal && *pending_request == request_id => {
                if durable_permit.is_some() {
                    return Err(SessionError::Persistence {
                        what: format!(
                            "submission {} already carries a durable claim permit",
                            request_id.content_hash
                        ),
                    });
                }
                *durable_permit = Some(permit);
                Ok(())
            }
            _ => Err(SessionError::Persistence {
                what: format!(
                    "submission {} lost its owned Pending state after durable claim commit",
                    request_id.content_hash
                ),
            }),
        }
    }

    /// Exactly-once execution under one typed request authority. Admission and
    /// causal meter commit have distinct ordinals; terminal publication and
    /// meter mutation occur in the same lock-held transition.
    ///
    /// # Errors
    /// Foreign/stale/conflicting authority, draining/paused gate state,
    /// capacity exhaustion, and corrupt session state fail closed. A panic or
    /// invalid returned charge becomes one replayable [`SubmitOutcome::Failed`]
    /// terminal rather than escaping as a partial mutation.
    #[allow(clippy::too_many_lines)]
    pub fn submit_once(
        &self,
        request_id: SubmissionRequestId,
        work: impl FnOnce() -> Charge,
    ) -> Result<SubmitOutcome, SessionError> {
        self.submit_once_inner(request_id, None, work)
    }

    /// Exactly-once execution with a durable pre-execution Pending claim.
    ///
    /// Only the call that inserts a fresh claim may invoke `work`. A recovered
    /// identical Pending claim returns [`SessionError::IndeterminateMutation`]
    /// because external side effects cannot be inferred. An existing terminal
    /// is recovered and replayed without invoking `work`. The terminal receipt
    /// and its audit event become restart-replayable when the scope is flushed.
    /// A new authority cannot run until the restarted governor has recovered
    /// every claim in its construction snapshot, across all sessions/scopes.
    ///
    /// # Errors
    /// A non-durable governor, foreign ledger, altered canonical program,
    /// Pending claim, claim conflict, admission refusal, execution failure, or
    /// corrupt durable state fails closed. Panics in `work` become a terminal
    /// [`SubmitOutcome::Failed`] as in [`Self::submit_once`].
    #[allow(clippy::too_many_arguments)]
    pub fn submit_once_durable(
        &self,
        ledger: &fs_ledger::Ledger,
        request_id: SubmissionRequestId,
        canonical_program: &str,
        work: impl FnOnce() -> Charge,
    ) -> Result<SubmitOutcome, SessionError> {
        let ledger_instance_id = self.recovery_ledger(ledger)?;
        let supplied_request_hash = bounded_request_digest(
            "idempotency_program_text_bytes",
            IDEMPOTENCY_PROGRAM_DOMAIN,
            canonical_program,
        )?;
        if supplied_request_hash != request_id.request_hash {
            return Err(SessionError::MutationConflict {
                kind: recovery::KIND_SUBMISSION,
                id: request_id.content_hash,
            });
        }
        let (open_id, token, gate, historical_generation) = {
            let inner = self.inner.lock().expect("governor lock");
            let open_id = inner
                .session_open_ids
                .get(&request_id.session.0)
                .copied()
                .ok_or(SessionError::UnknownSession {
                    id: request_id.session.0,
                })?;
            let token = inner.tokens.get(&request_id.session.0).cloned().ok_or(
                SessionError::UnknownSession {
                    id: request_id.session.0,
                },
            )?;
            let current_open = Self::current_open_identity(&inner, request_id.session)?;
            if current_open != request_id.session_open {
                return Err(SessionError::MutationAuthorityMismatch {
                    kind: recovery::KIND_SUBMISSION,
                    id: request_id.content_hash,
                });
            }
            let current_generation = inner
                .gate_generations
                .get(&request_id.session.0)
                .copied()
                .unwrap_or(0);
            if request_id.generation > current_generation {
                return Err(SessionError::StaleMutationGeneration {
                    kind: recovery::KIND_SUBMISSION,
                    id: request_id.session.0,
                    supplied: request_id.generation,
                    current: current_generation,
                });
            }
            (
                open_id,
                token,
                inner
                    .open_requests
                    .get(&open_id)
                    .and_then(|replay| replay.gate.clone()),
                request_id.generation < current_generation,
            )
        };
        self.recover_open(ledger, open_id, token.as_ref().clone(), gate)?;
        let claim_exists = ledger
            .session_mutation_claim(&request_id.content_hash)
            .map_err(|error| SessionError::Persistence {
                what: format!(
                    "durable submission claim preflight {} failed: {error}",
                    request_id.content_hash
                ),
            })?
            .is_some();
        if historical_generation || claim_exists {
            return self.recover_submission(ledger, request_id, canonical_program);
        }
        {
            let inner = self.inner.lock().expect("governor lock");
            self.ensure_durable_recovery_complete(&inner)?;
        }
        self.submit_once_inner(
            request_id,
            Some((ledger, ledger_instance_id, canonical_program)),
            work,
        )
    }

    #[allow(clippy::too_many_lines)]
    fn submit_once_inner(
        &self,
        request_id: SubmissionRequestId,
        durable: Option<(&fs_ledger::Ledger, fs_ledger::LedgerInstanceId, &str)>,
        work: impl FnOnce() -> Charge,
    ) -> Result<SubmitOutcome, SessionError> {
        if request_id.governor_id != self.id {
            return Err(SessionError::MutationAuthorityMismatch {
                kind: "submission-request",
                id: request_id.content_hash,
            });
        }
        let session = request_id.session;
        let (admission_ordinal, ledger_scope) = {
            let mut g = self.inner.lock().expect("governor lock");
            match g.idempotency.get(&request_id) {
                Some(IdemState::Done {
                    admission_ordinal,
                    receipt,
                    meter_receipt,
                    ..
                }) => {
                    return Ok(SubmitOutcome::Duplicate {
                        admission_ordinal: *admission_ordinal,
                        receipt: *receipt,
                        enforcement: meter_receipt.enforcement.clone(),
                        meter_receipt: meter_receipt.clone(),
                    });
                }
                Some(IdemState::Failed {
                    admission_ordinal,
                    receipt,
                    evidence,
                    ..
                }) => {
                    return Ok(SubmitOutcome::Failed {
                        admission_ordinal: *admission_ordinal,
                        receipt: *receipt,
                        evidence: evidence.as_ref().clone(),
                    });
                }
                Some(IdemState::Pending { .. }) => return Ok(SubmitOutcome::InFlight),
                None => {}
            }
            if self.durable_sink.is_some() && durable.is_none() {
                return Err(SessionError::DurableLedgerRequired {
                    kind: recovery::KIND_SUBMISSION,
                    authority: request_id.content_hash,
                });
            }
            let token = g
                .tokens
                .get(&session.0)
                .cloned()
                .ok_or(SessionError::UnknownSession { id: session.0 })?;
            let current_open = Self::current_open_identity(&g, session)?;
            if request_id.session_open != current_open {
                return Err(SessionError::MutationAuthorityMismatch {
                    kind: "submission-request",
                    id: request_id.content_hash,
                });
            }
            let current_generation = g.gate_generations.get(&session.0).copied().unwrap_or(0);
            if request_id.generation != current_generation {
                return Err(SessionError::StaleMutationGeneration {
                    kind: "submission-request",
                    id: session.0,
                    supplied: request_id.generation,
                    current: current_generation,
                });
            }
            if let Some(pending) = g.pending_pause.get(&session.0) {
                return Err(SessionError::PauseAlreadyPending {
                    id: session.0,
                    requested_ordinal: pending.request_id.requested_ordinal,
                });
            }
            if g.gate_phases.get(&session.0) == Some(&GatePhase::ReadyToResume) {
                return Err(SessionError::ResumeNotActivated {
                    id: session.0,
                    generation: current_generation,
                });
            }
            if g.gates
                .get(&session.0)
                .is_some_and(|gate| gate.is_requested())
            {
                return Err(SessionError::SessionGateDraining {
                    id: session.0,
                    generation: current_generation,
                });
            }
            let key_index = g
                .idempotency_keys
                .get(&session.0)
                .expect("open session has submission key index");
            if let Some(existing) = key_index.get(&request_id.key_hash)
                && existing != &request_id
            {
                return Err(SessionError::MutationConflict {
                    kind: "submission-request",
                    id: existing.content_hash,
                });
            }
            if key_index.len() >= MAX_IDEMPOTENCY_KEYS_PER_SESSION {
                return Err(SessionError::LimitExceeded {
                    resource: "idempotency_keys_per_session",
                    limit: MAX_IDEMPOTENCY_KEYS_PER_SESSION,
                    observed_at_least: key_index.len().saturating_add(1),
                });
            }
            let report_count = g.meter_report_ids.get(&session.0).map_or(0, BTreeSet::len);
            let reserved_reports = g
                .reserved_meter_reports
                .get(&session.0)
                .copied()
                .unwrap_or(0);
            if report_count.saturating_add(reserved_reports) >= MAX_METER_REPORTS_PER_SESSION {
                return Err(SessionError::LimitExceeded {
                    resource: "meter_reports_per_session",
                    limit: MAX_METER_REPORTS_PER_SESSION,
                    observed_at_least: report_count
                        .saturating_add(reserved_reports)
                        .saturating_add(1),
                });
            }
            let future_ordinals = g
                .reserved_meter_ordinals
                .checked_add(1)
                .and_then(|advance| u64::try_from(advance).ok())
                .and_then(|advance| g.next_meter_commit_ordinal.checked_add(advance))
                .ok_or(SessionError::LimitExceeded {
                    resource: "meter_commit_ordinal",
                    limit: i64::MAX as usize,
                    observed_at_least: usize::MAX,
                })?;
            if future_ordinals > i64::MAX as u64 {
                return Err(SessionError::LimitExceeded {
                    resource: "meter_commit_ordinal",
                    limit: i64::MAX as usize,
                    observed_at_least: i64::MAX as usize,
                });
            }
            let retained_bytes = SUBMISSION_REQUEST_RETAINED_BYTES
                .checked_add(MAX_IDEMPOTENCY_TERMINAL_RETAINED_BYTES)
                .and_then(|bytes| bytes.checked_add(MAX_METER_RECEIPT_RETAINED_BYTES))
                .ok_or(SessionError::LimitExceeded {
                    resource: "retained_bytes_per_scope",
                    limit: MAX_RETAINED_BYTES_PER_SCOPE,
                    observed_at_least: usize::MAX,
                })?;
            ensure_retained_capacity(&g, &token.ledger_scope, retained_bytes)?;
            let admission_ordinal =
                g.next_submission_ordinal
                    .checked_add(1)
                    .ok_or(SessionError::LimitExceeded {
                        resource: "submission_ordinal",
                        limit: i64::MAX as usize,
                        observed_at_least: usize::MAX,
                    })?;
            if admission_ordinal > i64::MAX as u64 {
                return Err(SessionError::LimitExceeded {
                    resource: "submission_ordinal",
                    limit: i64::MAX as usize,
                    observed_at_least: i64::MAX as usize,
                });
            }
            if let Some(existing) = g.submission_admissions.get(&admission_ordinal) {
                return Err(SessionError::Persistence {
                    what: format!(
                        "submission admission ordinal {admission_ordinal} is already owned by {}",
                        existing.content_hash
                    ),
                });
            }
            g.next_submission_ordinal = admission_ordinal;
            g.submission_admissions
                .insert(admission_ordinal, request_id);
            g.reserved_meter_ordinals += 1;
            *g.reserved_meter_reports
                .get_mut(&session.0)
                .expect("open session has report reservation count") += 1;
            *g.pending_submissions
                .get_mut(&session.0)
                .expect("open session has pending-submission count") += 1;
            g.idempotency.insert(
                request_id,
                IdemState::Pending {
                    admission_ordinal,
                    request_id,
                    reserved_terminal_bytes: MAX_IDEMPOTENCY_TERMINAL_RETAINED_BYTES,
                    reserved_meter_bytes: MAX_METER_RECEIPT_RETAINED_BYTES,
                    durable_permit: None,
                },
            );
            g.idempotency_keys
                .get_mut(&session.0)
                .expect("open session has submission key index")
                .insert(request_id.key_hash, request_id);
            commit_retained_bytes(&mut g, &token.ledger_scope, retained_bytes);
            (admission_ordinal, token.ledger_scope.clone())
        };

        if let Some((ledger, ledger_instance_id, canonical_program)) = durable {
            let payload = recovery::encode_submission_payload(request_id);
            let claim = fs_ledger::session_registry::SessionMutationClaim {
                authority: request_id.content_hash,
                ledger_instance_id,
                governor_hash: self.id,
                session_open_hash: request_id.session_open,
                kind: recovery::KIND_SUBMISSION,
                session: request_id.session.0,
                ledger_scope: &ledger_scope,
                generation: request_id.generation,
                causal_ordinal: Some(admission_ordinal),
                payload: &payload,
            };
            let claim_result = match ledger.claim_session_mutation(&claim) {
                Ok(result) => result,
                Err(error) => {
                    self.rollback_submission_admission(
                        request_id,
                        admission_ordinal,
                        &ledger_scope,
                    )?;
                    return match &error {
                        fs_ledger::LedgerError::Invalid { field, .. }
                            if field == "session_claim.authority" =>
                        {
                            Err(SessionError::MutationConflict {
                                kind: recovery::KIND_SUBMISSION,
                                id: request_id.content_hash,
                            })
                        }
                        _ => Err(SessionError::Persistence {
                            what: format!(
                                "durable submission claim {} failed before work: {error}",
                                request_id.content_hash
                            ),
                        }),
                    };
                }
            };
            match claim_result {
                fs_ledger::session_registry::SessionMutationClaimResult::Claimed { permit } => {
                    self.attach_submission_permit(request_id, admission_ordinal, permit)?;
                }
                fs_ledger::session_registry::SessionMutationClaimResult::Pending { .. } => {
                    self.rollback_submission_admission(
                        request_id,
                        admission_ordinal,
                        &ledger_scope,
                    )?;
                    return Err(SessionError::IndeterminateMutation {
                        kind: recovery::KIND_SUBMISSION,
                        authority: request_id.content_hash,
                    });
                }
                fs_ledger::session_registry::SessionMutationClaimResult::Terminal { .. } => {
                    self.rollback_submission_admission(
                        request_id,
                        admission_ordinal,
                        &ledger_scope,
                    )?;
                    return self.recover_submission(ledger, request_id, canonical_program);
                }
            }
        }

        let work_result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(work));
        let mut g = self.inner.lock().expect("governor lock");
        let (reserved_terminal_bytes, reserved_meter_bytes, durable_permit) = match g
            .idempotency
            .get(&request_id)
        {
            Some(IdemState::Pending {
                admission_ordinal: pending_ordinal,
                request_id: pending_request,
                reserved_terminal_bytes,
                reserved_meter_bytes,
                durable_permit,
            }) if *pending_ordinal == admission_ordinal && *pending_request == request_id => (
                *reserved_terminal_bytes,
                *reserved_meter_bytes,
                *durable_permit,
            ),
            Some(IdemState::Pending { .. } | IdemState::Done { .. } | IdemState::Failed { .. })
            | None => {
                return Err(SessionError::Persistence {
                    what: format!(
                        "session {} submission request lost its owned pending generation before terminal publication",
                        session.0
                    ),
                });
            }
        };

        let completion = match work_result {
            Ok(charge) => {
                let report_id = Self::submission_meter_report_id(request_id);
                match self.commit_meter_locked(&mut g, report_id, charge, true) {
                    Ok(meter_receipt) => SubmissionCompletion::Done(charge, meter_receipt),
                    Err(error) => {
                        *g.reserved_meter_reports
                            .get_mut(&session.0)
                            .expect("open session has report reservation count") -= 1;
                        g.reserved_meter_ordinals -= 1;
                        SubmissionCompletion::Failed(RetainedEvidence::capture(&error.to_string()))
                    }
                }
            }
            Err(payload) => {
                *g.reserved_meter_reports
                    .get_mut(&session.0)
                    .expect("open session has report reservation count") -= 1;
                g.reserved_meter_ordinals -= 1;
                SubmissionCompletion::Failed(panic_evidence(payload.as_ref()))
            }
        };
        let receipt = submission_receipt(request_id, &ledger_scope, admission_ordinal, &completion);
        let terminal_event_ordinal = match &completion {
            SubmissionCompletion::Done(_, meter_receipt) => meter_receipt.commit_ordinal,
            SubmissionCompletion::Failed(_) => admission_ordinal,
        };
        *g.pending_submissions
            .get_mut(&session.0)
            .expect("open session has pending-submission count") -= 1;
        let (outcome, terminal_retained_bytes, release_meter_bytes, terminal_succeeded) =
            match completion {
                SubmissionCompletion::Done(charge, meter_receipt) => {
                    let terminal_retained_bytes =
                        enforcement_retained_bytes(&meter_receipt.enforcement);
                    let enforcement = meter_receipt.enforcement.clone();
                    g.idempotency.insert(
                        request_id,
                        IdemState::Done {
                            admission_ordinal,
                            receipt,
                            charge,
                            meter_receipt: meter_receipt.clone(),
                            durable_permit,
                        },
                    );
                    (
                        SubmitOutcome::Executed {
                            admission_ordinal,
                            charge,
                            enforcement,
                            meter_receipt,
                            receipt,
                        },
                        terminal_retained_bytes,
                        0,
                        true,
                    )
                }
                SubmissionCompletion::Failed(evidence) => {
                    let terminal_retained_bytes = evidence.preview.len();
                    g.idempotency.insert(
                        request_id,
                        IdemState::Failed {
                            admission_ordinal,
                            receipt,
                            evidence: Arc::new(evidence.clone()),
                            durable_permit,
                        },
                    );
                    (
                        SubmitOutcome::Failed {
                            admission_ordinal,
                            receipt,
                            evidence,
                        },
                        terminal_retained_bytes,
                        reserved_meter_bytes,
                        false,
                    )
                }
            };
        if terminal_retained_bytes > reserved_terminal_bytes {
            return Err(SessionError::Persistence {
                what: format!(
                    "session {} terminal state requires {terminal_retained_bytes} retained bytes but reserved only {reserved_terminal_bytes}",
                    session.0
                ),
            });
        }
        bump_scope_revision(&mut g, &ledger_scope);
        let scope = g
            .scopes
            .get_mut(&ledger_scope)
            .expect("registered session scope");
        if terminal_succeeded {
            scope.dirty_causal.insert((
                terminal_event_ordinal,
                DirtyCausalMutation::Submission(request_id),
            ));
        } else {
            scope
                .dirty_submission_failures
                .insert((terminal_event_ordinal, request_id));
        }
        release_retained_bytes(
            &mut g,
            &ledger_scope,
            reserved_terminal_bytes - terminal_retained_bytes + release_meter_bytes,
        );
        Ok(outcome)
    }

    /// The canonical idempotency key: separately domain-hashed agent/program
    /// inputs plus their exact lengths under a final domain. Memory stays
    /// fixed-size after the bounded input hashes.
    ///
    /// # Errors
    /// [`SessionError::LimitExceeded`] before hashing when either input exceeds
    /// [`MAX_IDEMPOTENCY_INPUT_BYTES`].
    pub fn idempotency_key(agent_key: &str, program_text: &str) -> Result<String, SessionError> {
        for (resource, value) in [
            ("idempotency_agent_key_bytes", agent_key),
            ("idempotency_program_text_bytes", program_text),
        ] {
            if value.len() > MAX_IDEMPOTENCY_INPUT_BYTES {
                return Err(SessionError::LimitExceeded {
                    resource,
                    limit: MAX_IDEMPOTENCY_INPUT_BYTES,
                    observed_at_least: value.len(),
                });
            }
        }
        let agent_digest = fs_blake3::hash_domain(IDEMPOTENCY_AGENT_DOMAIN, agent_key.as_bytes());
        let program_digest =
            fs_blake3::hash_domain(IDEMPOTENCY_PROGRAM_DOMAIN, program_text.as_bytes());
        let mut payload = Vec::with_capacity(80);
        payload.extend_from_slice(
            &u64::try_from(agent_key.len())
                .expect("bounded idempotency agent key length fits u64")
                .to_le_bytes(),
        );
        payload.extend_from_slice(agent_digest.as_bytes());
        payload.extend_from_slice(
            &u64::try_from(program_text.len())
                .expect("bounded idempotency program length fits u64")
                .to_le_bytes(),
        );
        payload.extend_from_slice(program_digest.as_bytes());
        Ok(format!(
            "fs-session-idem-v3:{}",
            fs_blake3::hash_domain(IDEMPOTENCY_KEY_DOMAIN, &payload)
        ))
    }

    /// Mint a bounded authority for one declared pressure action in the
    /// session's current execution generation.
    ///
    /// # Errors
    /// Blank/oversized keys, unknown sessions, or corrupt immutable open state
    /// are refused without minting an authority.
    pub fn pressure_action_id(
        &self,
        session: SessionId,
        client_key: &str,
    ) -> Result<PressureActionId, SessionError> {
        let key_hash = bounded_request_digest(
            "pressure_action_key_bytes",
            PRESSURE_ACTION_ID_DOMAIN,
            client_key,
        )?;
        let (session_open, generation) = self.mutation_context(session)?;
        let mut payload = Vec::new();
        push_hash(&mut payload, self.id);
        payload.extend_from_slice(&session.0.to_le_bytes());
        push_hash(&mut payload, session_open);
        payload.extend_from_slice(&generation.to_le_bytes());
        push_hash(&mut payload, key_hash);
        Ok(PressureActionId {
            governor_id: self.id,
            session,
            session_open,
            generation,
            content_hash: fs_blake3::hash_domain(PRESSURE_ACTION_ID_DOMAIN, &payload),
        })
    }

    /// Apply memory pressure at `level` (1..=3 ONLY): ladder steps
    /// `1..=level` are emitted IN THE DECLARED ORDER with attribution. Spill
    /// and coarsen remain `Declared` orchestration work; this governor does not
    /// falsely mark their subsystem effects complete. The
    /// `PauseSerializeResume` step requests
    /// cancellation on the session's OWN gate, resolved by `SessionId`
    /// from the binding made at [`Governor::open_session_gated`] — no
    /// gate crosses this API, so pausing a different session's work is
    /// unrepresentable (bead gp3.13). The request event is phase
    /// `Requested`; it becomes `Complete` only through
    /// [`Governor::acknowledge_pause`] with a ledger-verified solver
    /// checkpoint receipt.
    ///
    /// # Errors
    /// - [`SessionError::InvalidPressureLevel`] for levels 0 and > 3.
    /// - [`SessionError::UnknownSession`].
    /// - [`SessionError::UngatedSession`] when level 3 targets a
    ///   session opened without a cancellation gate. Validation is
    ///   ATOMIC: no step fires and nothing is ledgered.
    /// - [`SessionError::PauseAlreadyPending`] while any earlier pause request
    ///   for the session remains unacknowledged.
    /// - [`SessionError::ResumeNotActivated`] while a fresh gate awaits explicit
    ///   resumed-worker activation.
    /// - [`SessionError::LimitExceeded`] for event or ordinal exhaustion.
    #[allow(clippy::too_many_lines)] // The ordered preflight and ladder commit are one state machine.
    pub fn apply_memory_pressure(
        &self,
        action_id: PressureActionId,
        level: u8,
    ) -> Result<PressureReceipt, SessionError> {
        if !(1..=3).contains(&level) {
            return Err(SessionError::InvalidPressureLevel { level });
        }
        if action_id.governor_id != self.id {
            return Err(SessionError::MutationAuthorityMismatch {
                kind: "pressure-action",
                id: action_id.content_hash,
            });
        }
        let session = action_id.session;
        let mut g = self.inner.lock().expect("governor lock");
        let ledger_scope = g
            .tokens
            .get(&session.0)
            .map(|token| token.ledger_scope.clone())
            .ok_or(SessionError::UnknownSession { id: session.0 })?;
        if let Some(replay) = g.pressure_actions.get(&action_id) {
            if replay.level != level {
                return Err(SessionError::MutationConflict {
                    kind: "pressure-action",
                    id: action_id.content_hash,
                });
            }
            let scope = g
                .scopes
                .get(&ledger_scope)
                .expect("registered session scope");
            let end = replay
                .event_start
                .checked_add(replay.event_len)
                .ok_or_else(|| SessionError::Persistence {
                    what: "pressure replay event range overflowed".to_string(),
                })?;
            let events = scope
                .events
                .get(replay.event_start..end)
                .ok_or_else(|| SessionError::Persistence {
                    what: "pressure replay event range is no longer retained".to_string(),
                })?
                .iter()
                .map(|event| event.as_ref().clone())
                .collect();
            return Ok(PressureReceipt {
                action_id,
                level,
                events,
                content_hash: replay.content_hash,
            });
        }
        self.ensure_durable_recovery_complete(&g)?;
        let current_open = Self::current_open_identity(&g, session)?;
        if action_id.session_open != current_open {
            return Err(SessionError::MutationAuthorityMismatch {
                kind: "pressure-action",
                id: action_id.content_hash,
            });
        }
        let current_generation = g.gate_generations.get(&session.0).copied().unwrap_or(0);
        if action_id.generation != current_generation {
            return Err(SessionError::StaleMutationGeneration {
                kind: "pressure-action",
                id: session.0,
                supplied: action_id.generation,
                current: current_generation,
            });
        }
        let action_count = g
            .pressure_action_ids
            .get(&session.0)
            .map_or(0, BTreeSet::len);
        if action_count >= MAX_PRESSURE_ACTIONS_PER_SESSION {
            return Err(SessionError::LimitExceeded {
                resource: "pressure_actions_per_session",
                limit: MAX_PRESSURE_ACTIONS_PER_SESSION,
                observed_at_least: action_count.saturating_add(1),
            });
        }
        if let Some(pending) = g.pending_pause.get(&session.0) {
            return Err(SessionError::PauseAlreadyPending {
                id: session.0,
                requested_ordinal: pending.request_id.requested_ordinal,
            });
        }
        if g.gate_phases.get(&session.0) == Some(&GatePhase::ReadyToResume) {
            let generation = *g
                .gate_generations
                .get(&session.0)
                .expect("ready gate has a generation");
            return Err(SessionError::ResumeNotActivated {
                id: session.0,
                generation,
            });
        }
        let scope = g
            .scopes
            .get(&ledger_scope)
            .expect("registered session scope");
        let reserve_completion = usize::from(usize::from(level) >= LADDER.len());
        let requested_event_count = scope
            .events
            .len()
            .saturating_add(scope.reserved_pause_completions)
            .saturating_add(usize::from(level))
            .saturating_add(reserve_completion);
        if requested_event_count > MAX_DEGRADATION_EVENTS_PER_SCOPE {
            return Err(SessionError::LimitExceeded {
                resource: "degradation_events_per_scope",
                limit: MAX_DEGRADATION_EVENTS_PER_SCOPE,
                observed_at_least: requested_event_count,
            });
        }
        let immediate_retained_bytes =
            LADDER[..usize::from(level)]
                .iter()
                .try_fold(0usize, |bytes, step| {
                    bytes
                        .checked_add(degradation_attribution(*step).len())
                        .ok_or(SessionError::LimitExceeded {
                            resource: "retained_bytes_per_scope",
                            limit: MAX_RETAINED_BYTES_PER_SCOPE,
                            observed_at_least: usize::MAX,
                        })
                })?;
        let reserved_completion_bytes = if usize::from(level) >= LADDER.len() {
            MAX_PAUSE_COMPLETION_RETAINED_BYTES
        } else {
            0
        };
        let retained_bytes = immediate_retained_bytes
            .checked_add(reserved_completion_bytes)
            .and_then(|bytes| bytes.checked_add(PRESSURE_ACTION_RETAINED_BYTES))
            .ok_or(SessionError::LimitExceeded {
                resource: "retained_bytes_per_scope",
                limit: MAX_RETAINED_BYTES_PER_SCOPE,
                observed_at_least: usize::MAX,
            })?;
        ensure_retained_capacity(&g, &ledger_scope, retained_bytes)?;
        let required_ordinal_advance = g
            .reserved_pause_ordinals
            .checked_add(usize::from(level))
            .and_then(|advance| advance.checked_add(reserve_completion))
            .and_then(|advance| i64::try_from(advance).ok())
            .ok_or(SessionError::LimitExceeded {
                resource: "degradation_ordinal",
                limit: i64::MAX as usize,
                observed_at_least: usize::MAX,
            })?;
        g.next_ordinal.checked_add(required_ordinal_advance).ok_or(
            SessionError::LimitExceeded {
                resource: "degradation_ordinal",
                limit: i64::MAX as usize,
                observed_at_least: usize::MAX,
            },
        )?;
        let final_ordinal = g
            .next_ordinal
            .checked_add(i64::from(level))
            .expect("reserved ordinal preflight covers immediate events");
        // Resolve the session's own gate BEFORE any step fires: a
        // refused level-3 request must not half-apply the ladder.
        let (gate, gate_generation) = if usize::from(level) >= LADDER.len() {
            let gate = g
                .gates
                .get(&session.0)
                .cloned()
                .ok_or(SessionError::UngatedSession { id: session.0 })?;
            let generation = *g
                .gate_generations
                .get(&session.0)
                .expect("gated session has a generation");
            generation
                .checked_add(1)
                .ok_or(SessionError::LimitExceeded {
                    resource: "pause_gate_generation",
                    limit: usize::MAX,
                    observed_at_least: usize::MAX,
                })?;
            if g.gate_phases.get(&session.0) != Some(&GatePhase::Running) {
                return Err(SessionError::Persistence {
                    what: format!(
                        "session {} has a gate and generation but no running gate phase",
                        session.0
                    ),
                });
            }
            (Some(gate), Some(generation))
        } else {
            (None, None)
        };
        if let Some(gate) = &gate {
            gate.request();
        }
        let first_ordinal = g.next_ordinal + 1;
        let pause_request_id = gate_generation.map(|generation| PauseRequestId {
            governor_id: self.id,
            session,
            gate_generation: generation,
            requested_ordinal: final_ordinal,
        });
        let mut fired = Vec::with_capacity(usize::from(level));
        for (i, step) in LADDER.iter().enumerate() {
            if i >= usize::from(level) {
                break;
            }
            let phase = match step {
                DegradationStep::SpillColdArenas | DegradationStep::CoarsenAdaptively => {
                    StepPhase::Declared
                }
                DegradationStep::PauseSerializeResume => StepPhase::Requested,
            };
            let is_pause = *step == DegradationStep::PauseSerializeResume;
            let event = DegradationEvent {
                session,
                step: *step,
                pressure_level: level,
                phase,
                attribution: degradation_attribution(*step).to_string(),
                ordinal: first_ordinal
                    + i64::try_from(i).expect("the fixed degradation ladder length fits i64"),
                requested_ordinal: None,
                checkpoint: None,
                gate_generation: if is_pause { gate_generation } else { None },
                pause_request_id: if is_pause { pause_request_id } else { None },
                pressure_action_id: Some(action_id),
            };
            fired.push(event.clone());
        }
        g.next_ordinal = final_ordinal;
        if let Some(requested) = fired
            .iter()
            .find(|event| event.phase == StepPhase::Requested)
        {
            g.pending_pause.insert(
                session.0,
                PendingPause {
                    request_id: requested
                        .pause_request_id
                        .expect("pause request carries acknowledgement authority"),
                    pressure_action_id: action_id,
                    reserved_retained_bytes: reserved_completion_bytes,
                },
            );
            g.reserved_pause_ordinals = g
                .reserved_pause_ordinals
                .checked_add(1)
                .expect("session bounds prevent pause-reservation overflow");
            g.scopes
                .get_mut(&ledger_scope)
                .expect("registered session scope")
                .reserved_pause_completions += 1;
        }
        let event_start = g
            .scopes
            .get(&ledger_scope)
            .expect("registered session scope")
            .events
            .len();
        {
            let scope = g
                .scopes
                .get_mut(&ledger_scope)
                .expect("registered session scope");
            scope.events.extend(fired.iter().cloned().map(Arc::new));
            scope
                .dirty_control
                .insert((first_ordinal, DirtyControlMutation::Pressure(action_id)));
        }
        let content_hash = pressure_receipt_hash(action_id, level, &fired);
        g.pressure_actions.insert(
            action_id,
            PressureReplay {
                level,
                event_start,
                event_len: fired.len(),
                content_hash,
            },
        );
        g.pressure_action_ids
            .get_mut(&session.0)
            .expect("open session has pressure-action index")
            .insert(action_id);
        commit_retained_bytes(&mut g, &ledger_scope, retained_bytes);
        bump_scope_revision(&mut g, &ledger_scope);
        Ok(PressureReceipt {
            action_id,
            level,
            events: fired,
            content_hash,
        })
    }

    /// Acknowledge one exact pending pause request with a ledger-minted solver
    /// checkpoint receipt: the ONLY route to a `Complete` pause event.
    /// The receipt is independently reverified against its physical ledger,
    /// persisted solver-state artifact, run provenance, and executor
    /// drain/finalize report before governor state is locked for completion.
    /// Identical replay of a completed request returns the same event and gate.
    /// If that gate was requested while it was still `ReadyToResume`, replay
    /// replaces the never-activated gate at the same generation and returns the
    /// replacement; the prior acknowledgement then fails exact-gate activation.
    /// Conflicting evidence or authority from another request fails closed.
    ///
    /// # Errors
    /// - [`SessionError::UnknownSession`].
    /// - [`SessionError::PauseRequestMismatch`] when the request is foreign,
    ///   stale, or does not name the session's pending/completed generation.
    /// - [`SessionError::PauseCheckpointMismatch`] for a forged, stale,
    ///   cross-session, cross-request, or unverified ledger receipt.
    /// - [`SessionError::PauseAcknowledgementConflict`] when a completed
    ///   request is replayed with different checkpoint evidence.
    /// - [`SessionError::LimitExceeded`] for event or ordinal exhaustion.
    #[allow(clippy::too_many_lines)] // One lock-held, rollback-free pause completion transition.
    pub fn acknowledge_pause(
        &self,
        request_id: PauseRequestId,
        ledger: &fs_ledger::Ledger,
        checkpoint: &fs_ledger::session_registry::SolverCheckpointReceipt,
    ) -> Result<PauseAcknowledgement, SessionError> {
        if request_id.governor_id != self.id {
            return Err(SessionError::PauseRequestMismatch {
                id: request_id.session.0,
                requested_ordinal: request_id.requested_ordinal,
            });
        }
        let session = request_id.session;
        let checkpoint_refusal = |reason| SessionError::PauseCheckpointMismatch {
            id: session.0,
            requested_ordinal: request_id.requested_ordinal,
            reason,
        };
        if checkpoint.session() != session.0 {
            return Err(checkpoint_refusal("cross-session"));
        }
        if checkpoint.gate_generation() != request_id.gate_generation {
            return Err(checkpoint_refusal("stale-generation"));
        }
        if checkpoint.pause_authority() != recovery::pause_ack_authority(request_id) {
            return Err(checkpoint_refusal("foreign-pause-request"));
        }
        ledger
            .verify_solver_checkpoint_receipt(checkpoint)
            .map_err(|_| checkpoint_refusal("unverified-ledger-receipt"))?;
        if self
            .durable_sink
            .is_some_and(|sink| sink != checkpoint.ledger_instance_id())
        {
            return Err(checkpoint_refusal("foreign-governor-ledger"));
        }
        let ledger_scope = {
            let g = self.inner.lock().expect("governor lock");
            let ledger_scope = g
                .tokens
                .get(&session.0)
                .map(|token| token.ledger_scope.clone())
                .ok_or(SessionError::UnknownSession { id: session.0 })?;
            if g.scopes
                .get(&ledger_scope)
                .and_then(|scope| scope.sink)
                .is_some_and(|sink| sink != checkpoint.ledger_instance_id())
            {
                return Err(checkpoint_refusal("foreign-scope-ledger"));
            }
            let pending = g
                .pending_pause
                .get(&session.0)
                .is_some_and(|pending| pending.request_id == request_id);
            let completed = g
                .completed_pause
                .get(&session.0)
                .is_some_and(|completed| completed.request_id == request_id);
            if !pending && !completed {
                return Err(SessionError::PauseRequestMismatch {
                    id: session.0,
                    requested_ordinal: request_id.requested_ordinal,
                });
            }
            ledger_scope
        };
        let evidence = RetainedEvidence::capture(&checkpoint.content_hash().to_hex());
        let mut g = self.inner.lock().expect("governor lock");
        if let Some(completed) = g
            .completed_pause
            .get(&session.0)
            .copied()
            .filter(|completed| completed.request_id == request_id)
        {
            if completed.checkpoint_byte_len != evidence.byte_len
                || completed.checkpoint_digest != evidence.digest
            {
                return Err(SessionError::PauseAcknowledgementConflict {
                    id: session.0,
                    requested_ordinal: request_id.requested_ordinal,
                });
            }
            let event = g
                .scopes
                .get(&ledger_scope)
                .expect("registered session scope")
                .events
                .get(completed.completion_event_index)
                .cloned()
                .ok_or_else(|| SessionError::Persistence {
                    what: format!(
                        "session {} completed pause event {} is missing from its scope",
                        session.0, completed.completion_ordinal
                    ),
                })?;
            if event.ordinal != completed.completion_ordinal
                || event.pause_request_id != Some(request_id)
            {
                return Err(SessionError::Persistence {
                    what: format!(
                        "session {} completed pause index no longer matches request ordinal {}",
                        session.0, request_id.requested_ordinal
                    ),
                });
            }
            let current_generation = g
                .gate_generations
                .get(&session.0)
                .copied()
                .ok_or(SessionError::UngatedSession { id: session.0 })?;
            if current_generation != completed.resume_generation {
                return Err(SessionError::PauseRequestMismatch {
                    id: session.0,
                    requested_ordinal: request_id.requested_ordinal,
                });
            }
            let mut resume_gate = g
                .gates
                .get(&session.0)
                .cloned()
                .ok_or(SessionError::UngatedSession { id: session.0 })?;
            if g.gate_phases.get(&session.0) == Some(&GatePhase::ReadyToResume)
                && resume_gate.is_requested()
            {
                // This gate was never activated, so it never named a running
                // execution generation. Replace only its Arc identity while
                // retaining the completion event and generation. The replayed
                // acknowledgement is then the sole activation authority; any
                // acknowledgement carrying the cancelled Arc fails ptr_eq.
                resume_gate = Arc::new(CancelGate::new());
                g.gates.insert(session.0, Arc::clone(&resume_gate));
            }
            return Ok(PauseAcknowledgement {
                request_id,
                event: event.as_ref().clone(),
                resume_gate,
                resume_generation: completed.resume_generation,
                gate_binding: completed.gate_binding,
                content_hash: completed.acknowledgement_hash,
            });
        }
        self.ensure_durable_recovery_complete(&g)?;
        let pending_submissions = g
            .pending_submissions
            .get(&session.0)
            .copied()
            .ok_or(SessionError::UnknownSession { id: session.0 })?;
        if pending_submissions != 0 {
            return Err(SessionError::PauseDrainPending {
                id: session.0,
                pending_submissions,
            });
        }
        let pending = *g
            .pending_pause
            .get(&session.0)
            .filter(|pending| pending.request_id == request_id)
            .ok_or(SessionError::PauseRequestMismatch {
                id: session.0,
                requested_ordinal: request_id.requested_ordinal,
            })?;
        let reserved = g
            .scopes
            .get(&ledger_scope)
            .expect("registered session scope")
            .reserved_pause_completions;
        if reserved == 0 || g.reserved_pause_ordinals == 0 {
            return Err(SessionError::Persistence {
                what: format!(
                    "session {} pending pause lost its completion row or ordinal reservation",
                    session.0
                ),
            });
        }
        let resume_generation =
            request_id
                .gate_generation
                .checked_add(1)
                .ok_or(SessionError::LimitExceeded {
                    resource: "pause_gate_generation",
                    limit: usize::MAX,
                    observed_at_least: usize::MAX,
                })?;
        let resume_gate = Arc::new(CancelGate::new());
        let event_count = g
            .scopes
            .get(&ledger_scope)
            .expect("registered session scope")
            .events
            .len();
        let reserved = g
            .scopes
            .get(&ledger_scope)
            .expect("registered session scope")
            .reserved_pause_completions;
        if reserved == 0 || event_count.saturating_add(reserved) > MAX_DEGRADATION_EVENTS_PER_SCOPE
        {
            return Err(SessionError::Persistence {
                what: format!("scope {ledger_scope:?} lost its reserved pause-completion capacity"),
            });
        }
        let ordinal = g
            .next_ordinal
            .checked_add(1)
            .ok_or(SessionError::LimitExceeded {
                resource: "degradation_ordinal",
                limit: i64::MAX as usize,
                observed_at_least: usize::MAX,
            })?;
        let event = DegradationEvent {
            session,
            step: DegradationStep::PauseSerializeResume,
            pressure_level: 3,
            phase: StepPhase::Complete,
            attribution: format!(
                "pause complete: ledger checkpoint {} binds solver-state artifact {}, run {}, \
                 and executor drain report {}; it acknowledges the request at ordinal {} and \
                 rotates gate generation {} to {resume_generation}",
                checkpoint.content_hash(),
                checkpoint.solver_state_artifact(),
                checkpoint.run().0,
                checkpoint.drain_report_hash(),
                request_id.requested_ordinal,
                request_id.gate_generation,
            ),
            ordinal,
            requested_ordinal: Some(request_id.requested_ordinal),
            checkpoint: Some(evidence),
            gate_generation: Some(request_id.gate_generation),
            pause_request_id: Some(request_id),
            pressure_action_id: Some(pending.pressure_action_id),
        };
        let event_retained_bytes = degradation_event_retained_bytes(&event)?;
        if event_retained_bytes > pending.reserved_retained_bytes {
            return Err(SessionError::Persistence {
                what: format!(
                    "pause completion requires {event_retained_bytes} retained bytes but generation {} reserved only {}",
                    request_id.gate_generation, pending.reserved_retained_bytes
                ),
            });
        }
        let released_retained_bytes = pending.reserved_retained_bytes - event_retained_bytes;
        let gate_binding = resumed_gate_binding(request_id, resume_generation);
        let acknowledgement_hash =
            pause_acknowledgement_hash(request_id, &event, resume_generation, gate_binding);
        g.pending_pause.remove(&session.0);
        g.reserved_pause_ordinals = g
            .reserved_pause_ordinals
            .checked_sub(1)
            .expect("pending pause owns one ordinal reservation");
        g.gates.insert(session.0, Arc::clone(&resume_gate));
        g.gate_generations.insert(session.0, resume_generation);
        g.gate_phases.insert(session.0, GatePhase::ReadyToResume);
        g.completed_pause.insert(
            session.0,
            CompletedPause {
                request_id,
                checkpoint_byte_len: event
                    .checkpoint
                    .as_ref()
                    .expect("completion carries evidence")
                    .byte_len,
                checkpoint_digest: event
                    .checkpoint
                    .as_ref()
                    .expect("completion carries evidence")
                    .digest,
                completion_event_index: event_count,
                completion_ordinal: ordinal,
                resume_generation,
                gate_binding,
                acknowledgement_hash,
            },
        );
        g.pause_acknowledgements.insert(
            request_id,
            PauseAcknowledgementReplay {
                completion_event_index: event_count,
                resume_generation,
                gate_binding,
                content_hash: acknowledgement_hash,
            },
        );
        g.next_ordinal = ordinal;
        {
            let scope = g
                .scopes
                .get_mut(&ledger_scope)
                .expect("registered session scope");
            scope.sink.get_or_insert(checkpoint.ledger_instance_id());
            scope.reserved_pause_completions -= 1;
            scope.events.push(Arc::new(event.clone()));
            scope.dirty_control.insert((
                ordinal,
                DirtyControlMutation::PauseAcknowledgement(request_id),
            ));
        }
        release_retained_bytes(&mut g, &ledger_scope, released_retained_bytes);
        bump_scope_revision(&mut g, &ledger_scope);
        Ok(PauseAcknowledgement {
            request_id,
            event,
            resume_gate,
            resume_generation,
            gate_binding,
            content_hash: acknowledgement_hash,
        })
    }

    /// Declare that resumed workers have adopted the acknowledgement's fresh
    /// gate. Identical activation is idempotent; pressure remains refused while
    /// the gate is only `ReadyToResume`.
    ///
    /// # Errors
    /// Foreign/stale acknowledgements, replaced gates, or a gate already
    /// requested before activation fail closed.
    pub fn activate_resume(
        &self,
        acknowledgement: &PauseAcknowledgement,
    ) -> Result<ResumeActivationReceipt, SessionError> {
        let request_id = acknowledgement.request_id;
        if request_id.governor_id != self.id {
            return Err(SessionError::ResumeAcknowledgementMismatch {
                id: request_id.session.0,
            });
        }
        let session = request_id.session.0;
        let mut g = self.inner.lock().expect("governor lock");
        if !g.tokens.contains_key(&session) {
            return Err(SessionError::UnknownSession { id: session });
        }
        let completed = g
            .completed_pause
            .get(&session)
            .copied()
            .ok_or(SessionError::ResumeAcknowledgementMismatch { id: session })?;
        let current_generation = g
            .gate_generations
            .get(&session)
            .copied()
            .ok_or(SessionError::UngatedSession { id: session })?;
        let current_gate = g
            .gates
            .get(&session)
            .cloned()
            .ok_or(SessionError::UngatedSession { id: session })?;
        let ledger_scope = g
            .tokens
            .get(&session)
            .expect("known session checked above")
            .ledger_scope
            .clone();
        let stored_event = g
            .scopes
            .get(&ledger_scope)
            .expect("registered session scope")
            .events
            .get(completed.completion_event_index)
            .ok_or_else(|| SessionError::Persistence {
                what: format!(
                    "session {session} completed pause event {} is missing from its scope",
                    completed.completion_ordinal
                ),
            })?;
        if stored_event.ordinal != completed.completion_ordinal
            || stored_event.pause_request_id != Some(completed.request_id)
        {
            return Err(SessionError::Persistence {
                what: format!(
                    "session {session} completed pause index no longer matches request ordinal {}",
                    completed.request_id.requested_ordinal
                ),
            });
        }
        if completed.request_id != request_id
            || stored_event.as_ref() != &acknowledgement.event
            || completed.resume_generation != acknowledgement.resume_generation
            || completed.gate_binding != acknowledgement.gate_binding
            || completed.acknowledgement_hash != acknowledgement.content_hash
            || current_generation != acknowledgement.resume_generation
            || !Arc::ptr_eq(&current_gate, &acknowledgement.resume_gate)
        {
            return Err(SessionError::ResumeAcknowledgementMismatch { id: session });
        }
        let session_open = Self::current_open_identity(&g, request_id.session)?;
        let activation_id = resume_activation_id(
            self.id,
            request_id.session,
            session_open,
            acknowledgement.content_hash,
            acknowledgement.resume_generation,
        );
        let receipt = resume_activation_receipt(
            activation_id,
            acknowledgement.content_hash,
            acknowledgement.gate_binding,
        );
        match g.gate_phases.get(&session).copied() {
            Some(GatePhase::ReadyToResume) => {
                self.ensure_durable_recovery_complete(&g)?;
                if current_gate.is_requested() {
                    return Err(SessionError::ResumeGateAlreadyRequested {
                        id: session,
                        generation: current_generation,
                    });
                }
                g.gate_phases.insert(session, GatePhase::Running);
                g.resume_activations.insert(activation_id, receipt);
                g.scopes
                    .get_mut(&ledger_scope)
                    .expect("registered session scope")
                    .dirty_control
                    .insert((
                        completed.completion_ordinal,
                        DirtyControlMutation::ResumeActivation(activation_id),
                    ));
                bump_scope_revision(&mut g, &ledger_scope);
                Ok(receipt)
            }
            Some(GatePhase::Running) => match g.resume_activations.get(&activation_id) {
                Some(stored) if *stored == receipt => Ok(receipt),
                Some(_) | None => Err(SessionError::Persistence {
                    what: format!(
                        "session {session} running generation {current_generation} lost its immutable resume-activation receipt"
                    ),
                }),
            },
            None => Err(SessionError::UngatedSession { id: session }),
        }
    }

    /// Whether a pause request is outstanding (requested, not yet
    /// acknowledged) for the session.
    ///
    /// # Errors
    /// [`SessionError::UnknownSession`].
    pub fn pause_pending(&self, session: SessionId) -> Result<bool, SessionError> {
        let g = self.inner.lock().expect("governor lock");
        if !g.tokens.contains_key(&session.0) {
            return Err(SessionError::UnknownSession { id: session.0 });
        }
        Ok(g.pending_pause.contains_key(&session.0))
    }

    /// Session consumption snapshot `(core_s, mem_peak, wall_s, throttled,
    /// paused)`.
    ///
    /// # Errors
    /// [`SessionError::UnknownSession`].
    pub fn consumption(
        &self,
        session: SessionId,
    ) -> Result<(f64, u64, f64, u32, u32), SessionError> {
        let g = self.inner.lock().expect("governor lock");
        let m = g
            .meters
            .get(&session.0)
            .ok_or(SessionError::UnknownSession { id: session.0 })?;
        Ok((m.core_s, m.mem_peak_bytes, m.wall_s, m.throttled, m.paused))
    }

    /// One bounded page of degradation events for the permit's exact scope.
    /// Results retain deterministic ordinal order.
    ///
    /// # Errors
    /// [`SessionError::ScopePermitMismatch`] for a permit minted by another
    /// governor, or [`SessionError::LimitExceeded`] when `limit` exceeds
    /// [`MAX_EVENT_PAGE_ROWS`].
    pub fn events_page(
        &self,
        permit: &ScopeFlushPermit,
        offset: usize,
        limit: usize,
    ) -> Result<Vec<DegradationEvent>, SessionError> {
        if permit.governor_id != self.id {
            return Err(SessionError::ScopePermitMismatch {
                scope: permit.ledger_scope.clone(),
            });
        }
        if limit > MAX_EVENT_PAGE_ROWS {
            return Err(SessionError::LimitExceeded {
                resource: "event_page_rows",
                limit: MAX_EVENT_PAGE_ROWS,
                observed_at_least: limit,
            });
        }
        let snapshot = {
            let g = self.inner.lock().expect("governor lock");
            let events = &g
                .scopes
                .get(&permit.ledger_scope)
                .ok_or_else(|| SessionError::UnknownLedgerScope {
                    scope: permit.ledger_scope.clone(),
                })?
                .events;
            let start = offset.min(events.len());
            let end = start.saturating_add(limit).min(events.len());
            events[start..end].to_vec()
        };
        self.wait_at_materialization_barrier_for_test();
        Ok(snapshot
            .into_iter()
            .map(|event| event.as_ref().clone())
            .collect())
    }

    /// Persist at most one bounded chunk for the permit's exact scope.
    /// Bounded source selection and cursor commit hold the governor mutex;
    /// evidence cloning, terminal encoding, and atomic database I/O do not.
    /// Call again while [`FlushReport::remaining_dirty`] is true.
    ///
    /// # Errors
    /// Foreign permits, concurrent same-scope flushes, sink mismatches,
    /// deterministic batch limits, explicit ledger transactions, and ledger
    /// failures are structured refusals. A failed append clears only the
    /// in-flight reservation and leaves every semantic cursor dirty.
    #[allow(clippy::too_many_lines)] // Explicit prepare / unlocked I/O / commit protocol.
    pub fn flush_scope_to_ledger(
        &self,
        permit: &ScopeFlushPermit,
        ledger: &fs_ledger::Ledger,
    ) -> Result<FlushReport, SessionError> {
        if permit.governor_id != self.id {
            return Err(SessionError::ScopePermitMismatch {
                scope: permit.ledger_scope.clone(),
            });
        }
        if ledger.in_transaction() {
            return Err(SessionError::Persistence {
                what: "session flush requires ownership of its atomic ledger transaction; an \
                       explicit transaction is already open and every flush cursor remains dirty"
                    .to_string(),
            });
        }
        let ledger_scope = permit.ledger_scope.clone();
        let sink_identity =
            ledger
                .checked_instance_id()
                .map_err(|error| SessionError::Persistence {
                    what: format!(
                        "ledger sink identity failed revalidation before scoped flush: {error}"
                    ),
                })?;
        if let Some(bound_sink) = self.durable_sink
            && bound_sink != sink_identity
        {
            return Err(SessionError::LedgerScopeSinkMismatch {
                scope: ledger_scope,
                bound_sink,
                attempted_sink: sink_identity,
            });
        }
        let snapshot = {
            let mut g = self.inner.lock().expect("governor lock");
            let scope =
                g.scopes
                    .get(&ledger_scope)
                    .ok_or_else(|| SessionError::UnknownLedgerScope {
                        scope: ledger_scope.clone(),
                    })?;
            if scope.in_flight.is_some() {
                return Err(SessionError::ScopeFlushInFlight {
                    scope: ledger_scope,
                });
            }
            if let Some(bound) = scope.sink
                && bound != sink_identity
            {
                return Err(SessionError::LedgerScopeSinkMismatch {
                    scope: ledger_scope,
                    bound_sink: bound,
                    attempted_sink: sink_identity,
                });
            }
            let generation =
                scope
                    .flush_generation
                    .checked_add(1)
                    .ok_or(SessionError::LimitExceeded {
                        resource: "scope_flush_generation",
                        limit: i64::MAX as usize,
                        observed_at_least: usize::MAX,
                    })?;
            let revision = scope.revision;
            // An open receipt is the durable authority prerequisite for every
            // later row belonging to that session.  Prioritize the bounded
            // open lane whenever it is dirty so lane rotation can never make
            // a dependent mutation visible in an earlier transaction.  If
            // the byte cap splits the open lane, the same override applies to
            // the next chunk; once it drains, rotation resumes at lane one.
            let start_flush_lane = if scope.dirty_open_receipts.is_empty() {
                scope.next_flush_lane
            } else {
                0
            };
            let next_flush_lane = (start_flush_lane + 1) % 4;
            let mut sources = Vec::with_capacity(MAX_FLUSH_ROWS.min(64));
            let mut selected_submissions = BTreeSet::new();

            // Rotate the first lane after every successful non-empty chunk.
            // This bounds preparation by dirty rows rather than retained state
            // and prevents continuously dirty report streams from starving
            // open, terminal, or degradation receipts.
            'lanes: for lane_offset in 0..4 {
                let remaining_rows = MAX_FLUSH_ROWS - sources.len();
                if remaining_rows == 0 {
                    break;
                }
                match (start_flush_lane + lane_offset) % 4 {
                    0 => {
                        let (dirty, has_more) = {
                            let scope = g.scopes.get(&ledger_scope).expect("scope checked above");
                            let dirty: Vec<SessionOpenId> = scope
                                .dirty_open_receipts
                                .iter()
                                .take(remaining_rows)
                                .copied()
                                .collect();
                            let has_more = scope.dirty_open_receipts.len() > dirty.len();
                            (dirty, has_more)
                        };
                        for open_id in dirty {
                            let replay = g.open_requests.get(&open_id).ok_or_else(|| {
                                SessionError::Persistence {
                                    what: format!(
                                        "scope dirty-open index references missing authority {}",
                                        open_id.content_hash
                                    ),
                                }
                            })?;
                            let receipt = Arc::clone(&replay.receipt);
                            let token = g.tokens.get(&open_id.session.0).cloned().ok_or_else(|| {
                                SessionError::Persistence {
                                    what: format!(
                                        "scope dirty-open index references missing token for session {}",
                                        open_id.session.0
                                    ),
                                }
                            })?;
                            sources.push(FlushSource::Open {
                                open_id,
                                receipt,
                                token,
                            });
                        }
                        if has_more {
                            break 'lanes;
                        }
                    }
                    1 => {
                        let (dirty, has_more) = {
                            let scope = g.scopes.get(&ledger_scope).expect("scope checked above");
                            let dirty: Vec<(u64, DirtyCausalMutation)> = scope
                                .dirty_causal
                                .iter()
                                .take(remaining_rows)
                                .copied()
                                .collect();
                            let has_more = scope.dirty_causal.len() > dirty.len();
                            (dirty, has_more)
                        };
                        for (indexed_ordinal, mutation) in dirty {
                            match mutation {
                                DirtyCausalMutation::Meter(report_id) => {
                                    let receipt =
                                        g.meter_reports.get(&report_id).ok_or_else(|| {
                                            SessionError::Persistence {
                                                what: format!(
                                                    "scope causal index references missing meter authority {}",
                                                    report_id.content_hash
                                                ),
                                            }
                                        })?;
                                    sources.push(FlushSource::Meter {
                                        indexed_ordinal,
                                        report_id,
                                        receipt: receipt.clone(),
                                    });
                                }
                                DirtyCausalMutation::Submission(request_id) => {
                                    let state =
                                        g.idempotency.get(&request_id).ok_or_else(|| {
                                            SessionError::Persistence {
                                                what: format!(
                                                    "scope causal index references missing submission {}",
                                                    request_id.content_hash
                                                ),
                                            }
                                        })?;
                                    sources.push(FlushSource::Submission {
                                        indexed_ordinal,
                                        request_id,
                                        state: state.clone(),
                                        failed_lane: false,
                                    });
                                    selected_submissions.insert(request_id);
                                }
                            }
                        }
                        if has_more {
                            break 'lanes;
                        }
                    }
                    2 => {
                        let (dirty, has_more) = {
                            let scope = g.scopes.get(&ledger_scope).expect("scope checked above");
                            let dirty: Vec<(u64, SubmissionRequestId)> = scope
                                .dirty_submission_failures
                                .iter()
                                .take(remaining_rows)
                                .copied()
                                .collect();
                            let has_more = scope.dirty_submission_failures.len() > dirty.len();
                            (dirty, has_more)
                        };
                        for (indexed_ordinal, request_id) in dirty {
                            let state = g.idempotency.get(&request_id).ok_or_else(|| {
                                SessionError::Persistence {
                                    what: format!(
                                        "scope failed-submission index references missing request {}",
                                        request_id.content_hash
                                    ),
                                }
                            })?;
                            sources.push(FlushSource::Submission {
                                indexed_ordinal,
                                request_id,
                                state: state.clone(),
                                failed_lane: true,
                            });
                            selected_submissions.insert(request_id);
                        }
                        if has_more {
                            break 'lanes;
                        }
                    }
                    3 => {
                        let (dirty, has_more) = {
                            let scope = g.scopes.get(&ledger_scope).expect("scope checked above");
                            let dirty: Vec<(i64, DirtyControlMutation)> = scope
                                .dirty_control
                                .iter()
                                .take(remaining_rows)
                                .copied()
                                .collect();
                            let has_more = scope.dirty_control.len() > dirty.len();
                            (dirty, has_more)
                        };
                        let mut deferred_pause = false;
                        for (indexed_ordinal, mutation) in dirty {
                            if let DirtyControlMutation::PauseAcknowledgement(request_id) = mutation
                            {
                                let scope =
                                    g.scopes.get(&ledger_scope).expect("scope checked above");
                                if has_dirty_submission_predecessor(
                                    scope,
                                    request_id,
                                    &selected_submissions,
                                ) {
                                    // The ledger rejects a pause acknowledgement while any
                                    // omitted draining-generation submission remains Pending.
                                    // A predecessor already selected into this atomic batch is
                                    // admissible; otherwise defer before sizing so lane rotation
                                    // cannot repeatedly prepare the same inadmissible prefix.
                                    deferred_pause = true;
                                    break;
                                }
                            }
                            let source = match mutation {
                                DirtyControlMutation::Pressure(action_id) => {
                                    let replay = g.pressure_actions.get(&action_id).ok_or_else(|| {
                                        SessionError::Persistence {
                                            what: format!(
                                                "scope control index references missing pressure action {}",
                                                action_id.content_hash
                                            ),
                                        }
                                    })?;
                                    let scope =
                                        g.scopes.get(&ledger_scope).expect("scope checked above");
                                    let events = scope
                                        .events
                                        .get(
                                            replay.event_start
                                                ..replay.event_start + replay.event_len,
                                        )
                                        .ok_or_else(|| SessionError::Persistence {
                                            what: format!(
                                                "pressure action {} lost its retained event group",
                                                action_id.content_hash
                                            ),
                                        })?
                                        .to_vec();
                                    FlushSource::Pressure {
                                        indexed_ordinal,
                                        action_id,
                                        replay: *replay,
                                        events,
                                    }
                                }
                                DirtyControlMutation::PauseAcknowledgement(request_id) => {
                                    let replay = g
                                        .pause_acknowledgements
                                        .get(&request_id)
                                        .ok_or_else(|| SessionError::Persistence {
                                            what: format!(
                                                "scope control index references missing pause acknowledgement {}",
                                                recovery::pause_ack_authority(request_id)
                                            ),
                                        })?;
                                    let event = g
                                        .scopes
                                        .get(&ledger_scope)
                                        .and_then(|scope| {
                                            scope.events.get(replay.completion_event_index)
                                        })
                                        .cloned()
                                        .ok_or_else(|| SessionError::Persistence {
                                            what: "pause acknowledgement lost its completion event"
                                                .to_string(),
                                        })?;
                                    let action_id = event.pressure_action_id.ok_or_else(|| {
                                        SessionError::Persistence {
                                            what: "pause completion lost its pressure action"
                                                .to_string(),
                                        }
                                    })?;
                                    let action_hash = g
                                        .pressure_actions
                                        .get(&action_id)
                                        .ok_or_else(|| SessionError::Persistence {
                                            what: "pause completion action receipt is missing"
                                                .to_string(),
                                        })?
                                        .content_hash;
                                    let gate = g.gates.get(&request_id.session.0).cloned().ok_or(
                                        SessionError::UngatedSession {
                                            id: request_id.session.0,
                                        },
                                    )?;
                                    FlushSource::PauseAcknowledgement {
                                        indexed_ordinal,
                                        request_id,
                                        replay: *replay,
                                        event,
                                        action_hash,
                                        gate,
                                        session_open: Self::current_open_identity(
                                            &g,
                                            request_id.session,
                                        )?,
                                    }
                                }
                                DirtyControlMutation::ResumeActivation(activation_id) => {
                                    let receipt = *g
                                        .resume_activations
                                        .get(&activation_id)
                                        .ok_or_else(|| SessionError::Persistence {
                                            what: format!(
                                                "scope control index references missing activation {}",
                                                activation_id.content_hash
                                            ),
                                        })?;
                                    FlushSource::ResumeActivation {
                                        indexed_ordinal,
                                        activation_id,
                                        receipt,
                                    }
                                }
                            };
                            sources.push(source);
                        }
                        if has_more && !deferred_pause {
                            break 'lanes;
                        }
                    }
                    _ => unreachable!("flush lane modulo four"),
                }
            }

            if sources.is_empty() {
                return Ok(FlushReport {
                    appended_rows: 0,
                    committed_terminals: 0,
                    encoded_bytes: 0,
                    remaining_dirty: false,
                });
            }
            let reservation_id =
                g.next_flush_reservation
                    .checked_add(1)
                    .ok_or(SessionError::LimitExceeded {
                        resource: "flush_reservation_ordinal",
                        limit: usize::MAX,
                        observed_at_least: usize::MAX,
                    })?;
            g.next_flush_reservation = reservation_id;
            g.scopes
                .get_mut(&ledger_scope)
                .expect("scope checked above")
                .in_flight = Some(reservation_id);
            FlushSnapshot {
                reservation_id,
                generation,
                revision,
                next_flush_lane,
                sources,
            }
        };

        self.wait_at_materialization_barrier_for_test();
        let snapshot_reservation_id = snapshot.reservation_id;
        let prepared = match snapshot.materialize(&ledger_scope) {
            Ok(prepared) => prepared,
            Err(error) => {
                let mut g = self.inner.lock().expect("governor lock");
                let scope = g
                    .scopes
                    .get_mut(&ledger_scope)
                    .expect("reserved scope remains registered");
                if scope.in_flight == Some(snapshot_reservation_id) {
                    scope.in_flight = None;
                }
                return Err(error);
            }
        };

        let event_groups: Vec<Vec<_>> = prepared
            .terminals
            .iter()
            .map(|terminal| {
                terminal
                    .events
                    .iter()
                    .map(BufferedLedgerEvent::as_row)
                    .collect()
            })
            .collect();
        let terminal_groups: Vec<_> = prepared
            .terminals
            .iter()
            .zip(&event_groups)
            .map(|(terminal, events)| {
                let claim = fs_ledger::session_registry::SessionMutationClaim {
                    authority: terminal.authority,
                    ledger_instance_id: sink_identity,
                    governor_hash: self.id,
                    session_open_hash: terminal.session_open,
                    kind: terminal.kind,
                    session: terminal.session.0,
                    ledger_scope: &ledger_scope,
                    generation: terminal.generation,
                    causal_ordinal: terminal.causal_ordinal,
                    payload: &terminal.payload,
                };
                fs_ledger::session_registry::SessionTerminalGroup {
                    terminal: fs_ledger::session_registry::SessionTerminalRow {
                        claim,
                        permit: terminal.permit,
                        receipt: &terminal.receipt,
                    },
                    events,
                }
            })
            .collect();
        let batch = fs_ledger::session_registry::SessionTerminalBatch {
            groups: &terminal_groups,
        };
        let persistence = ledger.append_session_terminal_batch(&batch);
        let (appended_rows, committed_terminals) = match persistence {
            Ok(fs_ledger::session_registry::SessionTerminalBatchResult::Committed {
                events_appended,
                ..
            }) => (events_appended, prepared.terminals.len()),
            Ok(fs_ledger::session_registry::SessionTerminalBatchResult::Replayed { .. }) => {
                (0, prepared.terminals.len())
            }
            Err(error) => {
                let mut g = self.inner.lock().expect("governor lock");
                let scope = g
                    .scopes
                    .get_mut(&ledger_scope)
                    .expect("reserved scope remains registered");
                if scope.in_flight == Some(prepared.reservation_id) {
                    scope.in_flight = None;
                }
                return Err(SessionError::Persistence {
                    what: format!(
                        "atomic bounded session batch failed; every semantic cursor remains dirty: {error}"
                    ),
                });
            }
        };

        let mut g = self.inner.lock().expect("governor lock");
        {
            let scope = g
                .scopes
                .get_mut(&ledger_scope)
                .expect("reserved scope remains registered");
            if scope.in_flight != Some(prepared.reservation_id) {
                return Err(SessionError::Persistence {
                    what: "scope flush reservation changed after a committed ledger batch; \
                           refusing to guess cursor ownership"
                        .to_string(),
                });
            }
            scope.in_flight = None;
            scope.sink.get_or_insert(sink_identity);
            scope.flush_generation = prepared.generation;
            scope.next_flush_lane = prepared.next_flush_lane;
        }
        for (open_id, content_hash) in prepared.open_marks {
            let still_current = g
                .open_requests
                .get(&open_id)
                .is_some_and(|replay| replay.receipt.content_hash == content_hash);
            if still_current {
                g.scopes
                    .get_mut(&ledger_scope)
                    .expect("committed scope remains registered")
                    .dirty_open_receipts
                    .remove(&open_id);
            }
        }
        for (report_id, content_hash) in prepared.meter_report_marks {
            let current = g.meter_reports.get(&report_id);
            if let Some(receipt) = current.filter(|receipt| receipt.content_hash == content_hash) {
                let commit_ordinal = receipt.commit_ordinal;
                g.scopes
                    .get_mut(&ledger_scope)
                    .expect("committed scope remains registered")
                    .dirty_causal
                    .remove(&(commit_ordinal, DirtyCausalMutation::Meter(report_id)));
            }
        }
        for (request_id, generation) in prepared.idempotency_marks {
            let still_current = match g.idempotency.get(&request_id) {
                Some(IdemState::Done {
                    admission_ordinal,
                    receipt,
                    meter_receipt,
                    ..
                }) => (*admission_ordinal, *receipt, meter_receipt.commit_ordinal) == generation,
                Some(IdemState::Failed {
                    admission_ordinal,
                    receipt,
                    ..
                }) => (*admission_ordinal, *receipt, *admission_ordinal) == generation,
                Some(IdemState::Pending { .. }) | None => false,
            };
            if still_current {
                let terminal_succeeded =
                    matches!(g.idempotency.get(&request_id), Some(IdemState::Done { .. }));
                let scope = g
                    .scopes
                    .get_mut(&ledger_scope)
                    .expect("committed scope remains registered");
                if terminal_succeeded {
                    scope
                        .dirty_causal
                        .remove(&(generation.2, DirtyCausalMutation::Submission(request_id)));
                } else {
                    scope
                        .dirty_submission_failures
                        .remove(&(generation.2, request_id));
                }
            }
        }
        for (ordinal, mutation, content_hash, owned_events) in prepared.control_marks {
            let still_current = match mutation {
                DirtyControlMutation::Pressure(action_id) => g
                    .pressure_actions
                    .get(&action_id)
                    .is_some_and(|replay| replay.content_hash == content_hash),
                DirtyControlMutation::PauseAcknowledgement(request_id) => g
                    .pause_acknowledgements
                    .get(&request_id)
                    .is_some_and(|replay| replay.content_hash == content_hash),
                DirtyControlMutation::ResumeActivation(activation_id) => g
                    .resume_activations
                    .get(&activation_id)
                    .is_some_and(|receipt| receipt.content_hash == content_hash),
            };
            if still_current {
                let scope = g
                    .scopes
                    .get_mut(&ledger_scope)
                    .expect("committed scope remains registered");
                scope.dirty_control.remove(&(ordinal, mutation));
                scope.flushed_events = scope.flushed_events.checked_add(owned_events).ok_or(
                    SessionError::LimitExceeded {
                        resource: "scope_event_cursor",
                        limit: usize::MAX,
                        observed_at_least: usize::MAX,
                    },
                )?;
                if scope.flushed_events > scope.events.len() {
                    return Err(SessionError::Persistence {
                        what: format!(
                            "scope committed-event count {} exceeds retained event count {}",
                            scope.flushed_events,
                            scope.events.len()
                        ),
                    });
                }
            }
        }
        let scope = g
            .scopes
            .get(&ledger_scope)
            .expect("committed scope remains registered");
        let remaining_dirty = scope.revision != prepared.revision
            || !scope.dirty_open_receipts.is_empty()
            || !scope.dirty_causal.is_empty()
            || !scope.dirty_submission_failures.is_empty()
            || !scope.dirty_control.is_empty();
        Ok(FlushReport {
            appended_rows,
            committed_terminals,
            encoded_bytes: prepared.encoded_bytes,
            remaining_dirty,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn identity_test_hash(byte: u8) -> fs_blake3::ContentHash {
        fs_blake3::hash_domain("fs-session-identity-test", &[byte])
    }

    fn solver_checkpoint(
        ledger: &fs_ledger::Ledger,
        request: PauseRequestId,
        label: &str,
    ) -> fs_ledger::session_registry::SolverCheckpointReceipt {
        let authority = request.checkpoint_authority();
        let mut run_bytes = [0_u8; 8];
        run_bytes.copy_from_slice(&authority.as_bytes()[..8]);
        let run = fs_exec::RunId(u64::from_le_bytes(run_bytes));
        let gate = fs_exec::CancelGate::new_clock_free();
        let tracker = fs_exec::DrainTracker::new(run, &gate);
        let worker = tracker.register_worker().expect("fixture worker");
        gate.request();
        drop(worker);
        let report = tracker.finalize().expect("fixture run drained");
        let snapshot =
            fs_exec::solver::envelope::seal(0x4653_434b_5054, 1, run.0, label.as_bytes());
        let artifact = ledger
            .put_artifact(
                fs_ledger::session_registry::SOLVER_STATE_ARTIFACT_KIND,
                &snapshot,
                None,
            )
            .expect("fixture solver-state artifact");
        ledger
            .attest_solver_checkpoint(fs_ledger::session_registry::SolverCheckpointClaim {
                session: request.session().0,
                pause_authority: authority,
                gate_generation: request.gate_generation(),
                solver_state_artifact: artifact.hash,
                drain_report: &report,
            })
            .expect("fixture checkpoint receipt")
    }

    #[test]
    fn governor_identity_versions_and_transports_fail_closed() {
        assert_eq!(DURABLE_GOVERNOR_IDENTITY_VERSION, 1);
        assert_eq!(SESSION_OPEN_IDENTITY_VERSION, 2);
        assert_eq!(SESSION_TOKEN_IDENTITY_VERSION, 1);
        assert_eq!(GATE_BINDING_IDENTITY_VERSION, 1);
        assert_eq!(METER_REPORT_IDENTITY_VERSION, 2);
        assert_eq!(PRESSURE_ACTION_IDENTITY_VERSION, 2);
        assert_eq!(SUBMISSION_AGENT_KEY_IDENTITY_VERSION, 1);
        assert_eq!(SUBMISSION_PROGRAM_IDENTITY_VERSION, 1);
        assert_eq!(SUBMISSION_REQUEST_IDENTITY_VERSION, 2);
        assert_eq!(RETAINED_EVIDENCE_IDENTITY_VERSION, 1);
        assert_eq!(PAUSE_ACKNOWLEDGEMENT_IDENTITY_VERSION, 1);
        assert_eq!(RESUME_ACTIVATION_IDENTITY_VERSION, 1);
        assert_eq!(SESSION_OPEN_RECEIPT_IDENTITY_VERSION, 1);
        assert_eq!(METER_RECEIPT_IDENTITY_VERSION, 1);
        assert_eq!(PRESSURE_RECEIPT_IDENTITY_VERSION, 1);
        assert_eq!(SUBMISSION_RECEIPT_IDENTITY_VERSION, 3);
        assert_eq!(PAUSE_ACKNOWLEDGEMENT_RECEIPT_IDENTITY_VERSION, 1);
        assert_eq!(RESUME_ACTIVATION_RECEIPT_IDENTITY_VERSION, 1);
        for (domain, version) in [
            (DURABLE_GOVERNOR_ID_DOMAIN, 1),
            (SESSION_OPEN_ID_DOMAIN, 2),
            (SESSION_TOKEN_IDENTITY_DOMAIN, 1),
            (GATE_BINDING_ID_DOMAIN, 1),
            (METER_REPORT_ID_DOMAIN, 2),
            (PRESSURE_ACTION_ID_DOMAIN, 2),
            (IDEMPOTENCY_AGENT_DOMAIN, 1),
            (IDEMPOTENCY_PROGRAM_DOMAIN, 1),
            (SUBMISSION_REQUEST_ID_DOMAIN, 2),
            (RETAINED_EVIDENCE_DOMAIN, 1),
            (PAUSE_ACKNOWLEDGEMENT_ID_DOMAIN, 1),
            (RESUME_ACTIVATION_ID_DOMAIN, 1),
            (SESSION_OPEN_RECEIPT_DOMAIN, 1),
            (METER_RECEIPT_DOMAIN, 1),
            (PRESSURE_RECEIPT_DOMAIN, 1),
            (SUBMISSION_RECEIPT_DOMAIN, 3),
            (PAUSE_ACKNOWLEDGEMENT_RECEIPT_DOMAIN, 1),
            (RESUME_ACTIVATION_RECEIPT_DOMAIN, 1),
        ] {
            assert!(domain.ends_with(&format!("v{version}")));
        }
        assert!(session_identity_transport_is_current(
            recovery::TERMINAL_SCHEMA_VERSION,
            None,
        ));
        for schema in [
            SESSION_OPEN_AUDIT_SCHEMA,
            SESSION_METER_AUDIT_SCHEMA,
            SESSION_SUBMISSION_AUDIT_SCHEMA,
            SESSION_DEGRADATION_AUDIT_SCHEMA,
        ] {
            assert!(session_identity_transport_is_current(
                recovery::TERMINAL_SCHEMA_VERSION,
                Some(schema),
            ));
        }
        assert!(!session_identity_transport_is_current(0, None));
        assert!(!session_identity_transport_is_current(
            recovery::TERMINAL_SCHEMA_VERSION + 1,
            None,
        ));
        assert!(!session_identity_transport_is_current(
            recovery::TERMINAL_SCHEMA_VERSION,
            Some("fs-session-open-v2"),
        ));
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn governor_authority_identity_fields_move_independently() {
        let h = identity_test_hash;

        let durable = DurableGovernorIdentitySource {
            ledger_instance_id: [1; 16],
            nonce: [2; 32],
        };
        let durable_hash = durable_governor_identity(&durable);
        let mut changed = durable.clone();
        changed.ledger_instance_id[0] ^= 1;
        assert_ne!(durable_hash, durable_governor_identity(&changed));
        changed = durable.clone();
        changed.nonce[0] ^= 1;
        assert_ne!(durable_hash, durable_governor_identity(&changed));

        let open = SessionOpenIdIdentitySource {
            governor_id: h(1),
            session: 2,
            client_key_digest: h(3),
        };
        let open_hash = session_open_authority_identity(&open);
        let mut changed = open.clone();
        changed.governor_id = h(4);
        assert_ne!(open_hash, session_open_authority_identity(&changed));
        changed = open.clone();
        changed.session += 1;
        assert_ne!(open_hash, session_open_authority_identity(&changed));
        changed = open.clone();
        changed.client_key_digest = h(5);
        assert_ne!(open_hash, session_open_authority_identity(&changed));

        let token = SessionTokenIdentitySource {
            session: 4,
            core_s_bits: 1.0_f64.to_bits(),
            mem_bytes: 5,
            wall_s_bits: 2.0_f64.to_bits(),
            cores: 6,
            ledger_scope: b"scope".to_vec(),
            operations: vec![b"a".to_vec(), b"b".to_vec()],
        };
        let token_hash = session_token_identity(&token);
        let mut changed = token.clone();
        changed.session += 1;
        assert_ne!(token_hash, session_token_identity(&changed));
        changed = token.clone();
        changed.core_s_bits ^= 1;
        assert_ne!(token_hash, session_token_identity(&changed));
        changed = token.clone();
        changed.mem_bytes += 1;
        assert_ne!(token_hash, session_token_identity(&changed));
        changed = token.clone();
        changed.wall_s_bits ^= 1;
        assert_ne!(token_hash, session_token_identity(&changed));
        changed = token.clone();
        changed.cores += 1;
        assert_ne!(token_hash, session_token_identity(&changed));
        changed = token.clone();
        changed.ledger_scope.push(b'x');
        assert_ne!(token_hash, session_token_identity(&changed));
        changed = token.clone();
        changed.operations.swap(0, 1);
        assert_ne!(token_hash, session_token_identity(&changed));
        changed = token.clone();
        changed.operations.push(b"c".to_vec());
        assert_ne!(token_hash, session_token_identity(&changed));

        let initial_gate = GateBindingIdentitySource::Session { open_id: h(7) };
        let initial_gate_hash = gate_binding_identity(&initial_gate);
        assert_ne!(
            initial_gate_hash,
            gate_binding_identity(&GateBindingIdentitySource::Session { open_id: h(8) })
        );
        let resumed_gate = GateBindingIdentitySource::Resumed {
            governor_id: h(9),
            session: 10,
            gate_generation: 11,
            requested_ordinal: 12,
            resume_generation: 13,
        };
        let resumed_gate_hash = gate_binding_identity(&resumed_gate);
        assert_ne!(initial_gate_hash, resumed_gate_hash);
        if let GateBindingIdentitySource::Resumed {
            governor_id,
            session,
            gate_generation,
            requested_ordinal,
            resume_generation,
        } = resumed_gate.clone()
        {
            for source in [
                GateBindingIdentitySource::Resumed {
                    governor_id: h(14),
                    session,
                    gate_generation,
                    requested_ordinal,
                    resume_generation,
                },
                GateBindingIdentitySource::Resumed {
                    governor_id,
                    session: session + 1,
                    gate_generation,
                    requested_ordinal,
                    resume_generation,
                },
                GateBindingIdentitySource::Resumed {
                    governor_id,
                    session,
                    gate_generation: gate_generation + 1,
                    requested_ordinal,
                    resume_generation,
                },
                GateBindingIdentitySource::Resumed {
                    governor_id,
                    session,
                    gate_generation,
                    requested_ordinal: requested_ordinal + 1,
                    resume_generation,
                },
                GateBindingIdentitySource::Resumed {
                    governor_id,
                    session,
                    gate_generation,
                    requested_ordinal,
                    resume_generation: resume_generation + 1,
                },
            ] {
                assert_ne!(resumed_gate_hash, gate_binding_identity(&source));
            }
        }

        let meter = MeterReportIdIdentitySource::Client {
            governor_id: h(15),
            session: 16,
            session_open: h(17),
            generation: 18,
            client_key_digest: h(19),
        };
        let meter_hash = meter_report_authority_identity(&meter);
        for source in [
            MeterReportIdIdentitySource::Client {
                governor_id: h(20),
                session: 16,
                session_open: h(17),
                generation: 18,
                client_key_digest: h(19),
            },
            MeterReportIdIdentitySource::Client {
                governor_id: h(15),
                session: 17,
                session_open: h(17),
                generation: 18,
                client_key_digest: h(19),
            },
            MeterReportIdIdentitySource::Client {
                governor_id: h(15),
                session: 16,
                session_open: h(18),
                generation: 18,
                client_key_digest: h(19),
            },
            MeterReportIdIdentitySource::Client {
                governor_id: h(15),
                session: 16,
                session_open: h(17),
                generation: 19,
                client_key_digest: h(19),
            },
            MeterReportIdIdentitySource::Client {
                governor_id: h(15),
                session: 16,
                session_open: h(17),
                generation: 18,
                client_key_digest: h(21),
            },
            MeterReportIdIdentitySource::Submission { request_id: h(22) },
        ] {
            assert_ne!(meter_hash, meter_report_authority_identity(&source));
        }

        let pressure = PressureActionIdIdentitySource {
            governor_id: h(23),
            session: 24,
            session_open: h(25),
            generation: 26,
            client_key_digest: h(27),
        };
        let pressure_hash = pressure_action_authority_identity(&pressure);
        let mut changed = pressure.clone();
        changed.governor_id = h(28);
        assert_ne!(pressure_hash, pressure_action_authority_identity(&changed));
        changed = pressure.clone();
        changed.session += 1;
        assert_ne!(pressure_hash, pressure_action_authority_identity(&changed));
        changed = pressure.clone();
        changed.session_open = h(29);
        assert_ne!(pressure_hash, pressure_action_authority_identity(&changed));
        changed = pressure.clone();
        changed.generation += 1;
        assert_ne!(pressure_hash, pressure_action_authority_identity(&changed));
        changed = pressure.clone();
        changed.client_key_digest = h(30);
        assert_ne!(pressure_hash, pressure_action_authority_identity(&changed));

        assert_ne!(
            submission_agent_key_identity(&SubmissionAgentKeyIdentitySource {
                value: b"agent-a".to_vec(),
            }),
            submission_agent_key_identity(&SubmissionAgentKeyIdentitySource {
                value: b"agent-b".to_vec(),
            }),
        );
        assert_ne!(
            submission_program_identity(&SubmissionProgramIdentitySource {
                value: b"program-a".to_vec(),
            }),
            submission_program_identity(&SubmissionProgramIdentitySource {
                value: b"program-b".to_vec(),
            }),
        );

        let submission = SubmissionRequestIdIdentitySource {
            governor_id: h(31),
            session: 32,
            session_open: h(33),
            generation: 34,
            key_hash: h(35),
            request_hash: h(36),
        };
        let submission_hash = submission_request_authority_identity(&submission);
        let mut changed = submission.clone();
        changed.governor_id = h(37);
        assert_ne!(
            submission_hash,
            submission_request_authority_identity(&changed)
        );
        changed = submission.clone();
        changed.session += 1;
        assert_ne!(
            submission_hash,
            submission_request_authority_identity(&changed)
        );
        changed = submission.clone();
        changed.session_open = h(38);
        assert_ne!(
            submission_hash,
            submission_request_authority_identity(&changed)
        );
        changed = submission.clone();
        changed.generation += 1;
        assert_ne!(
            submission_hash,
            submission_request_authority_identity(&changed)
        );
        changed = submission.clone();
        changed.key_hash = h(39);
        assert_ne!(
            submission_hash,
            submission_request_authority_identity(&changed)
        );

        assert_ne!(
            retained_evidence_identity(&RetainedEvidenceIdentitySource {
                value: b"evidence-a".to_vec(),
            }),
            retained_evidence_identity(&RetainedEvidenceIdentitySource {
                value: b"evidence-b".to_vec(),
            }),
        );

        let pause = PauseAcknowledgementIdIdentitySource {
            governor_id: h(40),
            session: 41,
            gate_generation: 42,
            requested_ordinal: 43,
        };
        let pause_hash = pause_acknowledgement_authority_identity(&pause);
        let mut changed = pause.clone();
        changed.governor_id = h(44);
        assert_ne!(
            pause_hash,
            pause_acknowledgement_authority_identity(&changed)
        );
        changed = pause.clone();
        changed.session += 1;
        assert_ne!(
            pause_hash,
            pause_acknowledgement_authority_identity(&changed)
        );
        changed = pause.clone();
        changed.gate_generation += 1;
        assert_ne!(
            pause_hash,
            pause_acknowledgement_authority_identity(&changed)
        );
        changed = pause.clone();
        changed.requested_ordinal += 1;
        assert_ne!(
            pause_hash,
            pause_acknowledgement_authority_identity(&changed)
        );

        let activation = ResumeActivationIdIdentitySource {
            governor_id: h(45),
            session: 46,
            session_open: h(47),
            acknowledgement_hash: h(48),
            resume_generation: 49,
        };
        let activation_hash = resume_activation_authority_identity(&activation);
        let mut changed = activation.clone();
        changed.governor_id = h(50);
        assert_ne!(
            activation_hash,
            resume_activation_authority_identity(&changed)
        );
        changed = activation.clone();
        changed.session += 1;
        assert_ne!(
            activation_hash,
            resume_activation_authority_identity(&changed)
        );
        changed = activation.clone();
        changed.session_open = h(51);
        assert_ne!(
            activation_hash,
            resume_activation_authority_identity(&changed)
        );
        changed = activation.clone();
        changed.acknowledgement_hash = h(52);
        assert_ne!(
            activation_hash,
            resume_activation_authority_identity(&changed)
        );
        changed = activation;
        changed.resume_generation += 1;
        assert_ne!(
            activation_hash,
            resume_activation_authority_identity(&changed)
        );

        assert_ne!(
            fs_blake3::hash_domain(SESSION_OPEN_ID_DOMAIN, b"canonical-order"),
            fs_blake3::hash_domain(SESSION_OPEN_ID_DOMAIN, b"order-canonical"),
        );
        assert_ne!(
            fs_blake3::hash_domain(SESSION_OPEN_ID_DOMAIN, b"domain"),
            fs_blake3::hash_domain(SESSION_TOKEN_IDENTITY_DOMAIN, b"domain"),
        );
    }

    fn receipt_test_evidence(byte: u8) -> RetainedEvidenceReceiptIdentitySource {
        RetainedEvidenceReceiptIdentitySource {
            preview: vec![byte],
            byte_len: 1,
            digest: identity_test_hash(byte),
        }
    }

    fn receipt_test_event(byte: u8) -> PressureEventIdentitySource {
        PressureEventIdentitySource {
            session: u64::from(byte),
            step_tag: 0,
            pressure_level: 1,
            phase_tag: 0,
            attribution: vec![byte],
            ordinal: i64::from(byte),
            requested_ordinal: Some(i64::from(byte) + 1),
            checkpoint: Some(receipt_test_evidence(byte)),
            gate_generation: Some(u64::from(byte) + 2),
            pause_request_id: None,
            pressure_action_id: None,
        }
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn governor_receipt_identity_fields_move_independently() {
        let h = identity_test_hash;

        let open = SessionOpenReceiptIdentitySource {
            open_id: h(1),
            token_digest: h(2),
            gate_identity: Some(h(3)),
            ledger_scope: b"scope".to_vec(),
        };
        let open_hash = session_open_receipt_identity(&open);
        let mut changed = open.clone();
        changed.open_id = h(4);
        assert_ne!(open_hash, session_open_receipt_identity(&changed));
        changed = open.clone();
        changed.token_digest = h(5);
        assert_ne!(open_hash, session_open_receipt_identity(&changed));
        changed = open.clone();
        changed.gate_identity = None;
        assert_ne!(open_hash, session_open_receipt_identity(&changed));
        changed = open.clone();
        changed.gate_identity = Some(h(6));
        assert_ne!(open_hash, session_open_receipt_identity(&changed));
        changed = open.clone();
        changed.ledger_scope.push(b'x');
        assert_ne!(open_hash, session_open_receipt_identity(&changed));

        let meter = MeterReceiptIdentitySource {
            report_id: h(7),
            commit_ordinal: 8,
            delta: ChargeIdentitySource {
                core_s_bits: 1.0_f64.to_bits(),
                mem_peak_bytes: 9,
                wall_s_bits: 2.0_f64.to_bits(),
            },
            before: MeterSnapshotIdentitySource {
                core_s_bits: 3.0_f64.to_bits(),
                mem_peak_bytes: 10,
                wall_s_bits: 4.0_f64.to_bits(),
                throttled: 11,
                paused: 12,
            },
            after: MeterSnapshotIdentitySource {
                core_s_bits: 5.0_f64.to_bits(),
                mem_peak_bytes: 13,
                wall_s_bits: 6.0_f64.to_bits(),
                throttled: 14,
                paused: 15,
            },
            enforcement: EnforcementIdentitySource::Ok,
        };
        let meter_hash = meter_receipt_identity(&meter);
        let mut changed = meter.clone();
        changed.report_id = h(16);
        assert_ne!(meter_hash, meter_receipt_identity(&changed));
        changed = meter.clone();
        changed.commit_ordinal += 1;
        assert_ne!(meter_hash, meter_receipt_identity(&changed));
        changed = meter.clone();
        changed.delta.core_s_bits ^= 1;
        assert_ne!(meter_hash, meter_receipt_identity(&changed));
        changed = meter.clone();
        changed.delta.mem_peak_bytes += 1;
        assert_ne!(meter_hash, meter_receipt_identity(&changed));
        changed = meter.clone();
        changed.delta.wall_s_bits ^= 1;
        assert_ne!(meter_hash, meter_receipt_identity(&changed));
        changed = meter.clone();
        changed.before.core_s_bits ^= 1;
        assert_ne!(meter_hash, meter_receipt_identity(&changed));
        changed = meter.clone();
        changed.before.mem_peak_bytes += 1;
        assert_ne!(meter_hash, meter_receipt_identity(&changed));
        changed = meter.clone();
        changed.before.wall_s_bits ^= 1;
        assert_ne!(meter_hash, meter_receipt_identity(&changed));
        changed = meter.clone();
        changed.before.throttled += 1;
        assert_ne!(meter_hash, meter_receipt_identity(&changed));
        changed = meter.clone();
        changed.before.paused += 1;
        assert_ne!(meter_hash, meter_receipt_identity(&changed));
        changed = meter.clone();
        changed.after.core_s_bits ^= 1;
        assert_ne!(meter_hash, meter_receipt_identity(&changed));
        changed = meter.clone();
        changed.after.mem_peak_bytes += 1;
        assert_ne!(meter_hash, meter_receipt_identity(&changed));
        changed = meter.clone();
        changed.after.wall_s_bits ^= 1;
        assert_ne!(meter_hash, meter_receipt_identity(&changed));
        changed = meter.clone();
        changed.after.throttled += 1;
        assert_ne!(meter_hash, meter_receipt_identity(&changed));
        changed = meter.clone();
        changed.after.paused += 1;
        assert_ne!(meter_hash, meter_receipt_identity(&changed));
        changed = meter.clone();
        changed.enforcement = EnforcementIdentitySource::Throttled {
            resource: b"core-seconds".to_vec(),
            used_bits: 7.0_f64.to_bits(),
            granted_bits: 8.0_f64.to_bits(),
        };
        let throttled_hash = meter_receipt_identity(&changed);
        assert_ne!(meter_hash, throttled_hash);
        let throttled = changed.clone();
        if let EnforcementIdentitySource::Throttled { resource, .. } = &mut changed.enforcement {
            resource.push(b'x');
        }
        assert_ne!(throttled_hash, meter_receipt_identity(&changed));
        changed = throttled.clone();
        if let EnforcementIdentitySource::Throttled { used_bits, .. } = &mut changed.enforcement {
            *used_bits ^= 1;
        }
        assert_ne!(throttled_hash, meter_receipt_identity(&changed));
        changed = throttled;
        if let EnforcementIdentitySource::Throttled { granted_bits, .. } = &mut changed.enforcement
        {
            *granted_bits ^= 1;
        }
        assert_ne!(throttled_hash, meter_receipt_identity(&changed));
        changed = meter.clone();
        changed.enforcement = EnforcementIdentitySource::Paused {
            resource: b"memory-bytes".to_vec(),
            used_bits: 9.0_f64.to_bits(),
            granted_bits: 10.0_f64.to_bits(),
            resume_hint: b"resume".to_vec(),
        };
        let paused_hash = meter_receipt_identity(&changed);
        if let EnforcementIdentitySource::Paused { resume_hint, .. } = &mut changed.enforcement {
            resume_hint.push(b'x');
        }
        assert_ne!(paused_hash, meter_receipt_identity(&changed));

        let mut event = receipt_test_event(17);
        let pressure = PressureReceiptIdentitySource {
            action_id: h(18),
            level: 2,
            events: vec![event.clone()],
        };
        let pressure_hash = pressure_receipt_identity(&pressure);
        let mut changed_pressure = pressure.clone();
        changed_pressure.action_id = h(19);
        assert_ne!(pressure_hash, pressure_receipt_identity(&changed_pressure));
        changed_pressure = pressure.clone();
        changed_pressure.level += 1;
        assert_ne!(pressure_hash, pressure_receipt_identity(&changed_pressure));
        changed_pressure = pressure.clone();
        changed_pressure.events.push(receipt_test_event(18));
        assert_ne!(pressure_hash, pressure_receipt_identity(&changed_pressure));
        changed_pressure.events.swap(0, 1);
        assert_ne!(pressure_hash, pressure_receipt_identity(&changed_pressure));
        let event_hash = pressure_receipt_identity(&pressure);
        event.session += 1;
        changed_pressure = pressure.clone();
        changed_pressure.events[0] = event.clone();
        assert_ne!(event_hash, pressure_receipt_identity(&changed_pressure));
        event = receipt_test_event(17);
        event.step_tag += 1;
        changed_pressure.events[0] = event.clone();
        assert_ne!(event_hash, pressure_receipt_identity(&changed_pressure));
        event = receipt_test_event(17);
        event.pressure_level += 1;
        changed_pressure.events[0] = event.clone();
        assert_ne!(event_hash, pressure_receipt_identity(&changed_pressure));
        event = receipt_test_event(17);
        event.phase_tag += 1;
        changed_pressure.events[0] = event.clone();
        assert_ne!(event_hash, pressure_receipt_identity(&changed_pressure));
        event = receipt_test_event(17);
        event.attribution.push(b'x');
        changed_pressure.events[0] = event.clone();
        assert_ne!(event_hash, pressure_receipt_identity(&changed_pressure));
        event = receipt_test_event(17);
        event.ordinal += 1;
        changed_pressure.events[0] = event.clone();
        assert_ne!(event_hash, pressure_receipt_identity(&changed_pressure));
        event = receipt_test_event(17);
        event.requested_ordinal = None;
        changed_pressure.events[0] = event.clone();
        assert_ne!(event_hash, pressure_receipt_identity(&changed_pressure));
        event = receipt_test_event(17);
        event.requested_ordinal = Some(99);
        changed_pressure.events[0] = event.clone();
        assert_ne!(event_hash, pressure_receipt_identity(&changed_pressure));
        event = receipt_test_event(17);
        event.checkpoint.as_mut().expect("checkpoint").byte_len += 1;
        changed_pressure.events[0] = event.clone();
        assert_ne!(event_hash, pressure_receipt_identity(&changed_pressure));
        event = receipt_test_event(17);
        event.checkpoint.as_mut().expect("checkpoint").digest = h(20);
        changed_pressure.events[0] = event.clone();
        assert_ne!(event_hash, pressure_receipt_identity(&changed_pressure));
        event = receipt_test_event(17);
        event.gate_generation = None;
        changed_pressure.events[0] = event.clone();
        assert_ne!(event_hash, pressure_receipt_identity(&changed_pressure));
        event = receipt_test_event(17);
        event.gate_generation = Some(100);
        changed_pressure.events[0] = event;
        assert_ne!(event_hash, pressure_receipt_identity(&changed_pressure));

        let done = SubmissionReceiptIdentitySource::Done {
            request_id: h(21),
            ledger_scope: b"scope".to_vec(),
            admission_ordinal: 22,
            charge: ChargeIdentitySource {
                core_s_bits: 1.0_f64.to_bits(),
                mem_peak_bytes: 23,
                wall_s_bits: 2.0_f64.to_bits(),
            },
            meter_receipt: h(24),
        };
        let done_hash = submission_receipt_identity(&done);
        let mut changed_done = done.clone();
        if let SubmissionReceiptIdentitySource::Done { request_id, .. } = &mut changed_done {
            *request_id = h(25);
        }
        assert_ne!(done_hash, submission_receipt_identity(&changed_done));
        changed_done = done.clone();
        if let SubmissionReceiptIdentitySource::Done { ledger_scope, .. } = &mut changed_done {
            ledger_scope.push(b'x');
        }
        assert_ne!(done_hash, submission_receipt_identity(&changed_done));
        changed_done = done.clone();
        if let SubmissionReceiptIdentitySource::Done {
            admission_ordinal, ..
        } = &mut changed_done
        {
            *admission_ordinal += 1;
        }
        assert_ne!(done_hash, submission_receipt_identity(&changed_done));
        changed_done = done.clone();
        if let SubmissionReceiptIdentitySource::Done { charge, .. } = &mut changed_done {
            charge.core_s_bits ^= 1;
        }
        assert_ne!(done_hash, submission_receipt_identity(&changed_done));
        changed_done = done.clone();
        if let SubmissionReceiptIdentitySource::Done { charge, .. } = &mut changed_done {
            charge.mem_peak_bytes += 1;
        }
        assert_ne!(done_hash, submission_receipt_identity(&changed_done));
        changed_done = done.clone();
        if let SubmissionReceiptIdentitySource::Done { charge, .. } = &mut changed_done {
            charge.wall_s_bits ^= 1;
        }
        assert_ne!(done_hash, submission_receipt_identity(&changed_done));
        changed_done = done.clone();
        if let SubmissionReceiptIdentitySource::Done { meter_receipt, .. } = &mut changed_done {
            *meter_receipt = h(26);
        }
        assert_ne!(done_hash, submission_receipt_identity(&changed_done));
        let failed = SubmissionReceiptIdentitySource::Failed {
            request_id: h(21),
            ledger_scope: b"scope".to_vec(),
            admission_ordinal: 22,
            evidence: receipt_test_evidence(27),
        };
        let failed_hash = submission_receipt_identity(&failed);
        assert_ne!(done_hash, failed_hash);
        let mut changed_failed = failed.clone();
        if let SubmissionReceiptIdentitySource::Failed { evidence, .. } = &mut changed_failed {
            evidence.byte_len += 1;
        }
        assert_ne!(failed_hash, submission_receipt_identity(&changed_failed));
        changed_failed = failed.clone();
        if let SubmissionReceiptIdentitySource::Failed { evidence, .. } = &mut changed_failed {
            evidence.digest = h(28);
        }
        assert_ne!(failed_hash, submission_receipt_identity(&changed_failed));

        let acknowledgement = PauseAcknowledgementReceiptIdentitySource {
            governor_id: h(29),
            session: 30,
            gate_generation: 31,
            requested_ordinal: 32,
            event: receipt_test_event(33),
            resume_generation: 34,
            gate_binding: h(35),
        };
        let acknowledgement_hash = pause_acknowledgement_receipt_identity(&acknowledgement);
        let mut changed = acknowledgement.clone();
        changed.governor_id = h(36);
        assert_ne!(
            acknowledgement_hash,
            pause_acknowledgement_receipt_identity(&changed)
        );
        changed = acknowledgement.clone();
        changed.session += 1;
        assert_ne!(
            acknowledgement_hash,
            pause_acknowledgement_receipt_identity(&changed)
        );
        changed = acknowledgement.clone();
        changed.gate_generation += 1;
        assert_ne!(
            acknowledgement_hash,
            pause_acknowledgement_receipt_identity(&changed)
        );
        changed = acknowledgement.clone();
        changed.requested_ordinal += 1;
        assert_ne!(
            acknowledgement_hash,
            pause_acknowledgement_receipt_identity(&changed)
        );
        changed = acknowledgement.clone();
        changed.event.ordinal += 1;
        assert_ne!(
            acknowledgement_hash,
            pause_acknowledgement_receipt_identity(&changed)
        );
        changed = acknowledgement.clone();
        changed.resume_generation += 1;
        assert_ne!(
            acknowledgement_hash,
            pause_acknowledgement_receipt_identity(&changed)
        );
        changed = acknowledgement.clone();
        changed.gate_binding = h(37);
        assert_ne!(
            acknowledgement_hash,
            pause_acknowledgement_receipt_identity(&changed)
        );

        let activation = ResumeActivationReceiptIdentitySource {
            activation_id: h(38),
            acknowledgement_hash: h(39),
            gate_binding: h(40),
        };
        let activation_hash = resume_activation_receipt_identity(&activation);
        let mut changed = activation.clone();
        changed.activation_id = h(41);
        assert_ne!(
            activation_hash,
            resume_activation_receipt_identity(&changed)
        );
        changed = activation.clone();
        changed.acknowledgement_hash = h(42);
        assert_ne!(
            activation_hash,
            resume_activation_receipt_identity(&changed)
        );
        changed = activation;
        changed.gate_binding = h(43);
        assert_ne!(
            activation_hash,
            resume_activation_receipt_identity(&changed)
        );

        assert_ne!(
            fs_blake3::hash_domain(SESSION_OPEN_RECEIPT_DOMAIN, b"canonical-order"),
            fs_blake3::hash_domain(SESSION_OPEN_RECEIPT_DOMAIN, b"order-canonical"),
        );
        assert_ne!(
            fs_blake3::hash_domain(SESSION_OPEN_RECEIPT_DOMAIN, b"domain"),
            fs_blake3::hash_domain(METER_RECEIPT_DOMAIN, b"domain"),
        );
    }

    #[test]
    fn governor_identity_nonsemantic_fields_do_not_move() {
        let h = identity_test_hash;
        let request = SubmissionRequestIdIdentitySource {
            governor_id: h(1),
            session: 2,
            session_open: h(3),
            generation: 4,
            key_hash: h(5),
            request_hash: h(6),
        };
        let mut changed_request = request.clone();
        changed_request.request_hash = h(7);
        assert_eq!(
            submission_request_authority_identity(&request),
            submission_request_authority_identity(&changed_request),
        );

        let event = receipt_test_event(8);
        let pressure = PressureReceiptIdentitySource {
            action_id: h(9),
            level: 1,
            events: vec![event.clone()],
        };
        let mut changed_pressure = pressure.clone();
        changed_pressure.events[0].pause_request_id = Some(PauseRequestId {
            governor_id: h(10),
            session: SessionId(11),
            gate_generation: 12,
            requested_ordinal: 13,
        });
        changed_pressure.events[0].pressure_action_id = Some(PressureActionId {
            governor_id: h(14),
            session: SessionId(15),
            session_open: h(16),
            generation: 17,
            content_hash: h(18),
        });
        changed_pressure.events[0]
            .checkpoint
            .as_mut()
            .expect("checkpoint")
            .preview
            .push(b'x');
        assert_eq!(
            pressure_receipt_identity(&pressure),
            pressure_receipt_identity(&changed_pressure),
        );

        let failure = SubmissionReceiptIdentitySource::Failed {
            request_id: h(19),
            ledger_scope: b"scope".to_vec(),
            admission_ordinal: 20,
            evidence: receipt_test_evidence(21),
        };
        let mut changed_failure = failure.clone();
        if let SubmissionReceiptIdentitySource::Failed { evidence, .. } = &mut changed_failure {
            evidence.preview.push(b'x');
        }
        assert_eq!(
            submission_receipt_identity(&failure),
            submission_receipt_identity(&changed_failure),
        );
    }

    struct LegacyGovernor {
        governor: super::Governor,
        next_mutation: AtomicU64,
    }

    use LegacyGovernor as Governor;

    impl LegacyGovernor {
        fn new() -> Self {
            Self {
                governor: super::Governor::new(),
                next_mutation: AtomicU64::new(1),
            }
        }

        fn next_key(&self, kind: &str, session: SessionId) -> String {
            let ordinal = self.next_mutation.fetch_add(1, Ordering::Relaxed);
            format!("legacy-test-{kind}-{}-{ordinal}", session.0)
        }

        fn open_session(&self, token: CapabilityToken) -> Result<ScopeFlushPermit, SessionError> {
            let open_id = self
                .governor
                .session_open_id(token.session, &self.next_key("open", token.session))?;
            self.governor
                .open_session(open_id, token)
                .map(|receipt| receipt.flush_permit())
        }

        fn open_session_gated(
            &self,
            token: CapabilityToken,
            gate: Arc<CancelGate>,
        ) -> Result<ScopeFlushPermit, SessionError> {
            let open_id = self
                .governor
                .session_open_id(token.session, &self.next_key("open", token.session))?;
            self.governor
                .open_session_gated(open_id, token, gate)
                .map(|receipt| receipt.flush_permit())
        }

        fn charge(&self, session: SessionId, delta: Charge) -> Result<Enforcement, SessionError> {
            let report_id = self
                .governor
                .meter_report_id(session, &self.next_key("meter", session))?;
            self.governor
                .charge(report_id, delta)
                .map(|receipt| receipt.enforcement.clone())
        }

        fn submit_once(
            &self,
            session: SessionId,
            key: &str,
            work: impl FnOnce() -> Charge,
        ) -> Result<SubmitOutcome, SessionError> {
            let request_id = self.governor.submission_request_id(session, key, key)?;
            self.governor.submit_once(request_id, work)
        }

        fn apply_memory_pressure(
            &self,
            session: SessionId,
            level: u8,
        ) -> Result<Vec<DegradationEvent>, SessionError> {
            let action_id = self
                .governor
                .pressure_action_id(session, &self.next_key("pressure", session))?;
            self.governor
                .apply_memory_pressure(action_id, level)
                .map(|receipt| receipt.events)
        }
    }

    impl core::ops::Deref for LegacyGovernor {
        type Target = super::Governor;

        fn deref(&self) -> &Self::Target {
            &self.governor
        }
    }

    fn test_token(session: u64, ledger_scope: &str) -> CapabilityToken {
        CapabilityToken {
            session: SessionId(session),
            ops: vec!["flux.*".to_string()],
            core_s: 1.0e9,
            mem_bytes: u64::MAX,
            wall_s: 1.0e9,
            cores: 1,
            ledger_scope: ledger_scope.to_string(),
        }
    }

    fn buffered(payload_len: usize) -> BufferedLedgerEvent {
        BufferedLedgerEvent {
            session: 1_u64.to_be_bytes(),
            t: 1,
            kind: "k",
            payload: "x".repeat(payload_len),
        }
    }

    fn durable_test_ledger_path(case: &str) -> String {
        static NEXT: AtomicU64 = AtomicU64::new(1);
        let ordinal = NEXT.fetch_add(1, Ordering::Relaxed);
        std::env::temp_dir()
            .join(format!(
                "fs-session-ujhp-{}-{ordinal}-{case}.ledger",
                std::process::id()
            ))
            .to_string_lossy()
            .into_owned()
    }

    fn bounded_terminal_with_events(event_count: usize) -> BufferedTerminal {
        let authority = fs_blake3::hash_domain(
            "org.frankensim.fs-session.test-terminal.v1",
            &u64::try_from(event_count)
                .expect("bounded test event count fits u64")
                .to_le_bytes(),
        );
        BufferedTerminal {
            authority,
            session_open: authority,
            kind: recovery::KIND_PRESSURE,
            session: SessionId(1),
            generation: 0,
            causal_ordinal: None,
            payload: Vec::new(),
            receipt: Vec::new(),
            events: (0..event_count)
                .map(|index| BufferedLedgerEvent {
                    session: 1_u64.to_be_bytes(),
                    t: i64::try_from(index).expect("bounded test ordinal fits i64"),
                    kind: "session.test-group",
                    payload: "{}".to_string(),
                })
                .collect(),
            permit: None,
        }
    }

    #[test]
    fn scoped_payload_preserves_and_escapes_the_exact_authority() {
        let payload = scoped_payload(
            "fs-session-test-v1",
            r#"alpha/"quoted"\branch"#,
            "\"value\":1",
        );
        assert_eq!(
            payload,
            r#"{"schema":"fs-session-test-v1","ledger_scope":"alpha/\"quoted\"\\branch","value":1}"#
        );
    }

    #[test]
    fn durable_governor_identity_reconstructs_only_on_the_same_ledger() {
        let ledger = fs_ledger::Ledger::open(":memory:").expect("fixture ledger");
        let nonce = DurableGovernorNonce::from_bytes([0x5a; 32]);
        let first = super::Governor::new_durable(&ledger, nonce).expect("first governor");
        let reopened =
            super::Governor::new_durable(&ledger, nonce).expect("reconstructed governor");
        assert_eq!(first.identity(), reopened.identity());
        let first_open = first
            .session_open_id(SessionId(30), "stable-open")
            .expect("first authority");
        let reopened_open = reopened
            .session_open_id(SessionId(30), "stable-open")
            .expect("reconstructed authority");
        assert_eq!(first_open, reopened_open);

        let foreign_ledger = fs_ledger::Ledger::open(":memory:").expect("foreign ledger");
        let foreign =
            super::Governor::new_durable(&foreign_ledger, nonce).expect("foreign durable governor");
        assert_ne!(first.identity(), foreign.identity());
        assert_ne!(
            first_open,
            foreign
                .session_open_id(SessionId(30), "stable-open")
                .expect("foreign authority")
        );
        assert_ne!(
            super::Governor::new().identity(),
            super::Governor::new().identity(),
            "ephemeral governor namespaces must not alias"
        );
    }

    #[test]
    fn submission_payload_carries_reconstructible_typed_authorities() {
        let governor = super::Governor::new();
        let session = SessionId(31);
        let open_id = governor
            .session_open_id(session, "payload-open")
            .expect("open authority");
        governor
            .open_session(open_id, test_token(session.0, "payload-authority"))
            .expect("fixture session");

        let success_id = governor
            .submission_request_id(session, "success-key", "success-program")
            .expect("success authority");
        governor
            .submit_once(success_id, Charge::default)
            .expect("success terminal");
        let failure_id = governor
            .submission_request_id(session, "failure-key", "failure-program")
            .expect("failure authority");
        governor
            .submit_once(failure_id, || Charge {
                core_s: f64::NAN,
                ..Charge::default()
            })
            .expect("invalid charge becomes a terminal failure");

        let inner = governor.inner.lock().expect("governor lock");
        let success_state = inner
            .idempotency
            .get(&success_id)
            .expect("success terminal retained");
        let (success, _) =
            buffered_submission_success("payload-authority", success_id, success_state)
                .expect("success payload");
        let derived_meter = super::Governor::submission_meter_report_id(success_id);
        for field in [
            format!("\"session_open\":\"{}\"", success_id.session_open),
            format!("\"generation\":{}", success_id.generation),
            format!("\"meter_report_id\":\"{}\"", derived_meter.content_hash),
        ] {
            assert!(
                success.payload.contains(&field),
                "success payload omitted {field}"
            );
        }

        let failure_state = inner
            .idempotency
            .get(&failure_id)
            .expect("failure terminal retained");
        let (failure, _) =
            buffered_submission_failure("payload-authority", failure_id, failure_state)
                .expect("failure payload");
        for field in [
            format!("\"session_open\":\"{}\"", failure_id.session_open),
            format!("\"generation\":{}", failure_id.generation),
        ] {
            assert!(
                failure.payload.contains(&field),
                "failure payload omitted {field}"
            );
        }
    }

    #[test]
    fn retained_evidence_bounds_preview_but_receipts_bind_the_full_tail() {
        let shared_prefix = "x".repeat(MAX_RETAINED_EVIDENCE_BYTES);
        let evidence_a = RetainedEvidence::capture(&format!("{shared_prefix}A"));
        let evidence_b = RetainedEvidence::capture(&format!("{shared_prefix}B"));
        assert_eq!(evidence_a.preview(), shared_prefix);
        assert_eq!(evidence_a.preview(), evidence_b.preview());
        assert_eq!(evidence_a.byte_len(), MAX_RETAINED_EVIDENCE_BYTES + 1);
        assert_ne!(evidence_a.digest(), evidence_b.digest());

        let governor = super::Governor::new();
        let open_id = governor
            .session_open_id(SessionId(1), "receipt-test-open")
            .expect("bounded open authority");
        let request_id = SubmissionRequestId {
            governor_id: governor.id,
            session: SessionId(1),
            session_open: open_id.content_hash,
            generation: 0,
            key_hash: fs_blake3::hash_domain(IDEMPOTENCY_AGENT_DOMAIN, b"key"),
            request_hash: fs_blake3::hash_domain(IDEMPOTENCY_PROGRAM_DOMAIN, b"program"),
            content_hash: fs_blake3::hash_domain(SUBMISSION_REQUEST_ID_DOMAIN, b"request"),
        };
        let receipt_a = submission_receipt(
            request_id,
            "scope",
            1,
            &SubmissionCompletion::Failed(evidence_a),
        );
        let receipt_b = submission_receipt(
            request_id,
            "scope",
            1,
            &SubmissionCompletion::Failed(evidence_b),
        );
        assert_ne!(
            receipt_a, receipt_b,
            "equal retained previews must not collapse distinct full evidence"
        );
    }

    #[test]
    fn registration_rebuilds_caller_controlled_spare_capacity() {
        let governor = Governor::new();
        let mut ledger_scope = String::with_capacity(4096);
        ledger_scope.push_str("canonical-capacity");
        let mut grant = String::with_capacity(4096);
        grant.push_str("flux.*");
        let mut ops = Vec::with_capacity(4096);
        ops.push(grant);
        let token = CapabilityToken {
            session: SessionId(88),
            ops,
            core_s: 1.0,
            mem_bytes: 1,
            wall_s: 1.0,
            cores: 1,
            ledger_scope,
        };
        governor.open_session(token).expect("valid bounded token");

        let inner = governor.inner.lock().expect("governor lock");
        let stored = &inner.tokens[&88];
        assert!(stored.ledger_scope.capacity() <= crate::token::MAX_LEDGER_SCOPE_BYTES);
        assert!(stored.ops.capacity() <= crate::token::MAX_CAPABILITY_OPS);
        assert!(stored.ops[0].capacity() <= crate::token::MAX_CAPABILITY_OP_BYTES);
    }

    #[test]
    fn retained_byte_budget_refuses_before_caller_work() {
        let governor = Governor::new();
        governor
            .open_session(test_token(89, "retained-budget"))
            .expect("fixture session");
        let key = "budget-key";
        let reservation = SUBMISSION_REQUEST_RETAINED_BYTES
            + MAX_IDEMPOTENCY_TERMINAL_RETAINED_BYTES
            + MAX_METER_RECEIPT_RETAINED_BYTES;
        {
            let mut inner = governor.inner.lock().expect("governor lock");
            inner.retained_bytes = MAX_RETAINED_BYTES_PER_SCOPE - reservation + 1;
            inner
                .scopes
                .get_mut("retained-budget")
                .expect("fixture scope")
                .retained_bytes = MAX_RETAINED_BYTES_PER_SCOPE - reservation + 1;
        }
        let executions = AtomicU64::new(0);
        assert!(matches!(
            governor.submit_once(SessionId(89), key, || {
                executions.fetch_add(1, Ordering::SeqCst);
                Charge::default()
            }),
            Err(SessionError::LimitExceeded {
                resource: "retained_bytes_per_scope",
                limit: MAX_RETAINED_BYTES_PER_SCOPE,
                ..
            })
        ));
        assert_eq!(executions.load(Ordering::SeqCst), 0);
        let inner = governor.inner.lock().expect("governor lock");
        assert!(inner.idempotency.is_empty());
        assert!(inner.idempotency_keys[&89].is_empty());
        assert_eq!(inner.next_submission_ordinal, 0);
    }

    #[test]
    fn pause_byte_budget_refuses_before_requesting_the_gate() {
        let governor = Governor::new();
        let gate = Arc::new(CancelGate::new());
        governor
            .open_session_gated(test_token(90, "pause-byte-budget"), Arc::clone(&gate))
            .expect("gated fixture session");
        {
            let mut inner = governor.inner.lock().expect("governor lock");
            inner.retained_bytes = MAX_RETAINED_BYTES_PER_SCOPE;
            inner
                .scopes
                .get_mut("pause-byte-budget")
                .expect("fixture scope")
                .retained_bytes = MAX_RETAINED_BYTES_PER_SCOPE;
        }
        assert!(matches!(
            governor.apply_memory_pressure(SessionId(90), 3),
            Err(SessionError::LimitExceeded {
                resource: "retained_bytes_per_scope",
                ..
            })
        ));
        assert!(!gate.is_requested());
        let inner = governor.inner.lock().expect("governor lock");
        assert!(inner.pending_pause.is_empty());
        assert_eq!(inner.reserved_pause_ordinals, 0);
        let scope = &inner.scopes["pause-byte-budget"];
        assert!(scope.events.is_empty());
        assert_eq!(scope.reserved_pause_completions, 0);
    }

    #[test]
    fn missing_gate_phase_refuses_before_requesting_the_gate() {
        let governor = Governor::new();
        let gate = Arc::new(CancelGate::new());
        governor
            .open_session_gated(test_token(91, "missing-gate-phase"), Arc::clone(&gate))
            .expect("gated fixture session");
        governor
            .inner
            .lock()
            .expect("governor lock")
            .gate_phases
            .remove(&91);

        assert!(matches!(
            governor.apply_memory_pressure(SessionId(91), 3),
            Err(SessionError::Persistence { what })
                if what.contains("no running gate phase")
        ));
        assert!(!gate.is_requested());
        let inner = governor.inner.lock().expect("governor lock");
        assert!(inner.pending_pause.is_empty());
        assert!(inner.scopes["missing-gate-phase"].events.is_empty());
    }

    #[test]
    fn externally_requested_registered_gate_enters_the_pause_protocol() {
        let governor = Governor::new();
        let gate = Arc::new(CancelGate::new());
        governor
            .open_session_gated(test_token(93, "external-cancel"), Arc::clone(&gate))
            .expect("gated fixture session");
        gate.request();

        let events = governor
            .apply_memory_pressure(SessionId(93), 3)
            .expect("an owned runtime cancellation can be formalized as a pause");
        let request_id = events
            .last()
            .and_then(|event| event.pause_request_id)
            .expect("pause request authority");
        let ledger = fs_ledger::Ledger::open(":memory:").expect("fixture ledger");
        let checkpoint = solver_checkpoint(&ledger, request_id, "external-cancel-checkpoint");
        let acknowledgement = governor
            .acknowledge_pause(request_id, &ledger, &checkpoint)
            .expect("pending cancellation can rotate to a fresh gate");
        governor
            .activate_resume(&acknowledgement)
            .expect("fresh gate activates");
        assert!(!acknowledgement.resume_gate().is_requested());
    }

    #[test]
    fn gate_generation_overflow_refuses_before_requesting_the_gate() {
        let governor = Governor::new();
        let gate = Arc::new(CancelGate::new());
        governor
            .open_session_gated(
                test_token(92, "gate-generation-overflow"),
                Arc::clone(&gate),
            )
            .expect("gated fixture session");
        governor
            .inner
            .lock()
            .expect("governor lock")
            .gate_generations
            .insert(92, u64::MAX);

        assert!(matches!(
            governor.apply_memory_pressure(SessionId(92), 3),
            Err(SessionError::LimitExceeded {
                resource: "pause_gate_generation",
                ..
            })
        ));
        assert!(!gate.is_requested());
        let inner = governor.inner.lock().expect("governor lock");
        assert!(inner.pending_pause.is_empty());
        assert_eq!(inner.reserved_pause_ordinals, 0);
        let scope = &inner.scopes["gate-generation-overflow"];
        assert!(scope.events.is_empty());
        assert_eq!(scope.reserved_pause_completions, 0);
    }

    #[test]
    fn bounded_flush_builder_enforces_exact_row_and_byte_limits() {
        let mut rows = Vec::new();
        let mut row_bytes = 0;
        for _ in 0..MAX_FLUSH_ROWS {
            assert!(push_bounded_event(&mut rows, &mut row_bytes, buffered(0)).unwrap());
        }
        let bytes_at_row_limit = row_bytes;
        assert!(!push_bounded_event(&mut rows, &mut row_bytes, buffered(0)).unwrap());
        assert_eq!(rows.len(), MAX_FLUSH_ROWS);
        assert_eq!(row_bytes, bytes_at_row_limit);

        let overhead = buffered(0).encoded_len().unwrap();
        let mut exact = Vec::new();
        let mut exact_bytes = 0;
        assert!(
            push_bounded_event(
                &mut exact,
                &mut exact_bytes,
                buffered(MAX_FLUSH_ENCODED_BYTES - overhead),
            )
            .unwrap()
        );
        assert_eq!(exact_bytes, MAX_FLUSH_ENCODED_BYTES);
        assert!(!push_bounded_event(&mut exact, &mut exact_bytes, buffered(0)).unwrap());
        assert_eq!(exact.len(), 1);
        assert!(matches!(
            push_bounded_event(
                &mut Vec::new(),
                &mut 0,
                buffered(MAX_FLUSH_ENCODED_BYTES - overhead + 1),
            ),
            Err(SessionError::LimitExceeded {
                resource: "flush_row_encoded_bytes",
                limit: MAX_FLUSH_ENCODED_BYTES,
                observed_at_least,
            }) if observed_at_least == MAX_FLUSH_ENCODED_BYTES + 1
        ));
    }

    #[test]
    fn event_and_ordinal_caps_refuse_before_mutation() {
        let governor = Governor::new();
        let permit = governor
            .open_session(test_token(1, "bounded"))
            .expect("fixture session");
        let fixture = DegradationEvent {
            session: SessionId(1),
            step: DegradationStep::SpillColdArenas,
            pressure_level: 1,
            phase: StepPhase::Declared,
            attribution: String::new(),
            ordinal: 0,
            requested_ordinal: None,
            checkpoint: None,
            gate_generation: None,
            pause_request_id: None,
            pressure_action_id: None,
        };
        {
            let mut inner = governor.inner.lock().expect("governor lock");
            inner
                .scopes
                .get_mut("bounded")
                .expect("fixture scope")
                .events = vec![Arc::new(fixture); MAX_DEGRADATION_EVENTS_PER_SCOPE - 1];
        }
        governor
            .apply_memory_pressure(SessionId(1), 1)
            .expect("exact event boundary is admitted");
        let before_refusal = {
            let inner = governor.inner.lock().expect("governor lock");
            (
                inner.next_ordinal,
                inner.scopes["bounded"].events.len(),
                inner.scopes["bounded"].revision,
            )
        };
        assert!(matches!(
            governor.apply_memory_pressure(SessionId(1), 1),
            Err(SessionError::LimitExceeded {
                resource: "degradation_events_per_scope",
                limit: MAX_DEGRADATION_EVENTS_PER_SCOPE,
                observed_at_least,
            }) if observed_at_least == MAX_DEGRADATION_EVENTS_PER_SCOPE + 1
        ));
        let after_refusal = {
            let inner = governor.inner.lock().expect("governor lock");
            (
                inner.next_ordinal,
                inner.scopes["bounded"].events.len(),
                inner.scopes["bounded"].revision,
            )
        };
        assert_eq!(after_refusal, before_refusal);
        assert_eq!(
            governor
                .events_page(&permit, MAX_DEGRADATION_EVENTS_PER_SCOPE - 1, 1)
                .expect("last event page")
                .len(),
            1
        );

        let ordinal_governor = Governor::new();
        ordinal_governor
            .open_session(test_token(2, "ordinal"))
            .expect("ordinal fixture session");
        {
            let mut inner = ordinal_governor.inner.lock().expect("governor lock");
            inner.next_ordinal = i64::MAX;
            inner.next_submission_ordinal = i64::MAX as u64;
        }
        assert!(matches!(
            ordinal_governor.apply_memory_pressure(SessionId(2), 1),
            Err(SessionError::LimitExceeded {
                resource: "degradation_ordinal",
                ..
            })
        ));
        let ran = AtomicU64::new(0);
        assert!(matches!(
            ordinal_governor.submit_once(SessionId(2), "ordinal-overflow", || {
                ran.fetch_add(1, Ordering::SeqCst);
                Charge::default()
            }),
            Err(SessionError::LimitExceeded {
                resource: "submission_ordinal",
                ..
            })
        ));
        assert_eq!(ran.load(Ordering::SeqCst), 0);
        let inner = ordinal_governor.inner.lock().expect("governor lock");
        assert!(inner.scopes["ordinal"].events.is_empty());
        assert!(inner.idempotency.is_empty());
        assert!(inner.idempotency_keys[&2].is_empty());
    }

    #[test]
    fn level_three_reserves_its_mandatory_completion_capacity() {
        let governor = Governor::new();
        let gate = Arc::new(CancelGate::new());
        governor
            .open_session_gated(test_token(9, "pause-capacity"), Arc::clone(&gate))
            .expect("gated fixture session");
        governor
            .open_session(test_token(10, "pause-capacity"))
            .expect("competing fixture session in the same scope");
        let fixture = DegradationEvent {
            session: SessionId(9),
            step: DegradationStep::SpillColdArenas,
            pressure_level: 1,
            phase: StepPhase::Declared,
            attribution: String::new(),
            ordinal: 0,
            requested_ordinal: None,
            checkpoint: None,
            gate_generation: None,
            pause_request_id: None,
            pressure_action_id: None,
        };
        {
            let mut inner = governor.inner.lock().expect("governor lock");
            inner
                .scopes
                .get_mut("pause-capacity")
                .expect("fixture scope")
                .events = vec![Arc::new(fixture); MAX_DEGRADATION_EVENTS_PER_SCOPE - 4];
        }
        let requested = governor
            .apply_memory_pressure(SessionId(9), 3)
            .expect("three events plus one completion reservation fit exactly");
        let request_id = requested
            .last()
            .and_then(|event| event.pause_request_id)
            .expect("level three mints request authority");
        {
            let inner = governor.inner.lock().expect("governor lock");
            let scope = &inner.scopes["pause-capacity"];
            assert_eq!(scope.events.len(), MAX_DEGRADATION_EVENTS_PER_SCOPE - 1);
            assert_eq!(scope.reserved_pause_completions, 1);
        }
        assert!(matches!(
            governor.apply_memory_pressure(SessionId(10), 1),
            Err(SessionError::LimitExceeded {
                resource: "degradation_events_per_scope",
                ..
            })
        ));
        let ledger = fs_ledger::Ledger::open(":memory:").expect("fixture ledger");
        let checkpoint = solver_checkpoint(&ledger, request_id, "checkpoint-at-cap");
        governor
            .acknowledge_pause(request_id, &ledger, &checkpoint)
            .expect("reserved completion remains admissible");
        let inner = governor.inner.lock().expect("governor lock");
        let scope = &inner.scopes["pause-capacity"];
        assert_eq!(scope.events.len(), MAX_DEGRADATION_EVENTS_PER_SCOPE);
        assert_eq!(scope.reserved_pause_completions, 0);
    }

    #[test]
    fn level_three_reserves_its_mandatory_completion_ordinal() {
        let governor = Governor::new();
        governor
            .open_session_gated(test_token(10, "pause-ordinal"), Arc::new(CancelGate::new()))
            .expect("gated fixture session");
        governor
            .open_session(test_token(11, "pause-ordinal"))
            .expect("interleaving fixture session");
        {
            let mut inner = governor.inner.lock().expect("governor lock");
            inner.next_ordinal = i64::MAX - 4;
        }

        let requested = governor
            .apply_memory_pressure(SessionId(10), 3)
            .expect("three immediate ordinals plus completion fit exactly");
        let request_id = requested
            .last()
            .and_then(|event| event.pause_request_id)
            .expect("level three mints request authority");
        {
            let inner = governor.inner.lock().expect("governor lock");
            assert_eq!(inner.next_ordinal, i64::MAX - 1);
            assert_eq!(inner.reserved_pause_ordinals, 1);
        }
        assert!(matches!(
            governor.apply_memory_pressure(SessionId(11), 1),
            Err(SessionError::LimitExceeded {
                resource: "degradation_ordinal",
                ..
            })
        ));
        let ledger = fs_ledger::Ledger::open(":memory:").expect("fixture ledger");
        let checkpoint = solver_checkpoint(&ledger, request_id, "checkpoint-at-ordinal-cap");
        governor
            .acknowledge_pause(request_id, &ledger, &checkpoint)
            .expect("reserved completion ordinal remains admissible");
        let inner = governor.inner.lock().expect("governor lock");
        assert_eq!(inner.next_ordinal, i64::MAX);
        assert_eq!(inner.reserved_pause_ordinals, 0);
    }

    #[test]
    fn identical_pause_acknowledgements_commit_once_and_replay() {
        let governor = Governor::new();
        governor
            .open_session_gated(test_token(12, "pause-replay"), Arc::new(CancelGate::new()))
            .expect("gated fixture session");
        let request_id = governor
            .apply_memory_pressure(SessionId(12), 3)
            .expect("pause request")
            .last()
            .and_then(|event| event.pause_request_id)
            .expect("request authority");
        let ledger = fs_ledger::Ledger::open(":memory:").expect("fixture ledger");
        let checkpoint = solver_checkpoint(&ledger, request_id, "same-checkpoint");
        let first = governor
            .acknowledge_pause(request_id, &ledger, &checkpoint)
            .expect("commit");
        let second = governor
            .acknowledge_pause(request_id, &ledger, &checkpoint)
            .expect("exact replay");
        assert_eq!(first.event, second.event);
        assert_eq!(first.resume_generation, second.resume_generation);
        assert!(Arc::ptr_eq(&first.resume_gate(), &second.resume_gate()));

        let inner = governor.inner.lock().expect("governor lock");
        let scope = &inner.scopes["pause-replay"];
        assert_eq!(
            scope
                .events
                .iter()
                .filter(|event| event.phase == StepPhase::Complete)
                .count(),
            1
        );
        assert_eq!(scope.reserved_pause_completions, 0);
        assert_eq!(inner.reserved_pause_ordinals, 0);
    }

    #[test]
    fn altered_pause_acknowledgement_cannot_activate_resume() {
        let governor = Governor::new();
        governor
            .open_session_gated(
                test_token(15, "pause-ack-integrity"),
                Arc::new(CancelGate::new()),
            )
            .expect("gated fixture session");
        let request_id = governor
            .apply_memory_pressure(SessionId(15), 3)
            .expect("pause request")
            .last()
            .and_then(|event| event.pause_request_id)
            .expect("request authority");
        let ledger = fs_ledger::Ledger::open(":memory:").expect("fixture ledger");
        let checkpoint = solver_checkpoint(&ledger, request_id, "checkpoint-claim");
        let acknowledgement = governor
            .acknowledge_pause(request_id, &ledger, &checkpoint)
            .expect("pause acknowledgement");
        let mut altered = acknowledgement.clone();
        altered.event.phase = StepPhase::Requested;

        assert_eq!(
            governor.activate_resume(&altered),
            Err(SessionError::ResumeAcknowledgementMismatch { id: 15 })
        );
        governor
            .activate_resume(&acknowledgement)
            .expect("unaltered acknowledgement remains authoritative");
    }

    #[test]
    fn independent_pending_pauses_share_and_consume_scope_reservations() {
        let governor = Governor::new();
        let mut requests = Vec::new();
        for session in [13, 14] {
            governor
                .open_session_gated(
                    test_token(session, "parallel-pauses"),
                    Arc::new(CancelGate::new()),
                )
                .expect("gated fixture session");
            requests.push(
                governor
                    .apply_memory_pressure(SessionId(session), 3)
                    .expect("parallel pause request")
                    .last()
                    .and_then(|event| event.pause_request_id)
                    .expect("request authority"),
            );
        }
        {
            let inner = governor.inner.lock().expect("governor lock");
            assert_eq!(inner.reserved_pause_ordinals, 2);
            assert_eq!(
                inner.scopes["parallel-pauses"].reserved_pause_completions,
                2
            );
        }
        let ledger = fs_ledger::Ledger::open(":memory:").expect("fixture ledger");
        for (index, request) in requests.into_iter().rev().enumerate() {
            let checkpoint = solver_checkpoint(&ledger, request, &format!("checkpoint-{index}"));
            governor
                .acknowledge_pause(request, &ledger, &checkpoint)
                .expect("each reservation is independently consumable");
        }
        let inner = governor.inner.lock().expect("governor lock");
        assert_eq!(inner.reserved_pause_ordinals, 0);
        assert_eq!(
            inner.scopes["parallel-pauses"].reserved_pause_completions,
            0
        );
    }

    #[test]
    fn same_scope_flush_reservation_refuses_a_race() {
        let governor = Governor::new();
        let permit = governor
            .open_session(test_token(3, "reserved"))
            .expect("fixture session");
        {
            let mut inner = governor.inner.lock().expect("governor lock");
            inner
                .scopes
                .get_mut("reserved")
                .expect("fixture scope")
                .in_flight = Some(7);
        }
        let ledger = fs_ledger::Ledger::open(":memory:").expect("fixture ledger");
        assert_eq!(
            governor.flush_scope_to_ledger(&permit, &ledger),
            Err(SessionError::ScopeFlushInFlight {
                scope: "reserved".to_string(),
            })
        );
        assert_eq!(ledger.table_count("events").unwrap(), 0);
    }

    #[test]
    #[allow(clippy::too_many_lines)] // One barrier protocol proves both bounded materialization surfaces.
    fn maximum_page_and_flush_materialize_without_holding_the_governor_lock() {
        let page_governor = Arc::new(Governor::new());
        let page_permit = page_governor
            .open_session(test_token(1, "materialize-page"))
            .expect("page fixture session");
        page_governor
            .open_session(test_token(2, "materialize-page-hot"))
            .expect("unrelated page hot-path session");
        let evidence = RetainedEvidence::capture(&"e".repeat(MAX_RETAINED_EVIDENCE_BYTES));
        let event = Arc::new(DegradationEvent {
            session: SessionId(1),
            step: DegradationStep::PauseSerializeResume,
            pressure_level: 3,
            phase: StepPhase::Complete,
            attribution: "maximum evidence-bearing page fixture".to_string(),
            ordinal: 1,
            requested_ordinal: Some(1),
            checkpoint: Some(evidence),
            gate_generation: Some(1),
            pause_request_id: None,
            pressure_action_id: None,
        });
        page_governor
            .inner
            .lock()
            .expect("governor lock")
            .scopes
            .get_mut("materialize-page")
            .expect("page scope")
            .events = vec![event; MAX_EVENT_PAGE_ROWS];

        let page_barrier = Arc::new(std::sync::Barrier::new(2));
        page_governor.set_materialization_barrier(Some(Arc::clone(&page_barrier)));
        let page_worker = {
            let governor = Arc::clone(&page_governor);
            std::thread::spawn(move || {
                governor
                    .events_page(&page_permit, 0, MAX_EVENT_PAGE_ROWS)
                    .expect("maximum page materializes")
            })
        };
        page_barrier.wait();
        let (page_hot_tx, page_hot_rx) = std::sync::mpsc::channel();
        let page_hot_worker = {
            let governor = Arc::clone(&page_governor);
            std::thread::spawn(move || {
                governor
                    .charge(SessionId(2), Charge::default())
                    .expect("unrelated charge while page materialization is paused");
                governor
                    .submit_once(SessionId(2), "page-hot-submit", Charge::default)
                    .expect("unrelated submission while page materialization is paused");
                page_hot_tx.send(()).expect("page hot-path witness");
            })
        };
        let page_hot = page_hot_rx.recv_timeout(std::time::Duration::from_secs(1));
        page_barrier.wait();
        page_governor.set_materialization_barrier(None);
        let page = page_worker.join().expect("page materialization worker");
        page_hot_worker.join().expect("page hot-path worker");
        assert!(
            page_hot.is_ok(),
            "unrelated scope charge/submit must finish while maximum page cloning is barrier-paused"
        );
        assert_eq!(page.len(), MAX_EVENT_PAGE_ROWS);

        let flush_governor = Arc::new(Governor::new());
        let mut flush_permit = None;
        for session in 0..MAX_SESSIONS_PER_SCOPE {
            let permit = flush_governor
                .open_session(test_token(
                    u64::try_from(session).expect("fixture id fits"),
                    "materialize-flush",
                ))
                .expect("maximum flush source row");
            flush_permit.get_or_insert(permit);
        }
        flush_governor
            .open_session(test_token(10_000, "materialize-flush-hot"))
            .expect("unrelated flush hot-path session");
        let flush_barrier = Arc::new(std::sync::Barrier::new(2));
        flush_governor.set_materialization_barrier(Some(Arc::clone(&flush_barrier)));
        let flush_worker = {
            let governor = Arc::clone(&flush_governor);
            let permit = flush_permit.expect("flush permit");
            std::thread::spawn(move || {
                let ledger = fs_ledger::Ledger::open(":memory:").expect("fixture ledger");
                governor
                    .flush_scope_to_ledger(&permit, &ledger)
                    .expect("maximum flush materializes")
            })
        };
        flush_barrier.wait();
        let (flush_hot_tx, flush_hot_rx) = std::sync::mpsc::channel();
        let flush_hot_worker = {
            let governor = Arc::clone(&flush_governor);
            std::thread::spawn(move || {
                governor
                    .charge(SessionId(10_000), Charge::default())
                    .expect("unrelated charge while flush materialization is paused");
                governor
                    .submit_once(SessionId(10_000), "flush-hot-submit", Charge::default)
                    .expect("unrelated submission while flush materialization is paused");
                flush_hot_tx.send(()).expect("flush hot-path witness");
            })
        };
        let flush_hot = flush_hot_rx.recv_timeout(std::time::Duration::from_secs(1));
        flush_barrier.wait();
        flush_governor.set_materialization_barrier(None);
        let report = flush_worker.join().expect("flush materialization worker");
        flush_hot_worker.join().expect("flush hot-path worker");
        assert!(
            flush_hot.is_ok(),
            "unrelated scope charge/submit must finish while maximum flush encoding is barrier-paused"
        );
        assert_eq!(report.appended_rows, MAX_FLUSH_ROWS);
        assert_eq!(report.committed_terminals, MAX_FLUSH_ROWS);
    }

    #[test]
    fn rotating_dirty_lanes_prevent_meter_starvation() {
        let governor = Governor::new();
        let mut permit = None;
        for session in 0..MAX_SESSIONS_PER_SCOPE {
            let opened = governor
                .open_session(test_token(
                    u64::try_from(session).expect("fixture id fits"),
                    "fair",
                ))
                .expect("fixture session");
            permit.get_or_insert(opened);
        }
        governor
            .submit_once(SessionId(0), "terminal", || Charge {
                core_s: 1.0,
                ..Charge::default()
            })
            .expect("terminal fixture");
        governor
            .apply_memory_pressure(SessionId(0), 1)
            .expect("event fixture");
        {
            let inner = governor.inner.lock().expect("governor lock");
            let scope = &inner.scopes["fair"];
            assert_eq!(
                scope.dirty_causal.len(),
                1,
                "the submission terminal occupies its meter-commit position"
            );
            assert!(scope.dirty_submission_failures.is_empty());
        }
        let ledger = fs_ledger::Ledger::open(":memory:").expect("fixture ledger");
        let permit = permit.expect("at least one fixture session");
        let first = governor
            .flush_scope_to_ledger(&permit, &ledger)
            .expect("open-receipt chunk");
        assert_eq!(first.appended_rows, MAX_FLUSH_ROWS);
        assert!(first.remaining_dirty);

        // Add one later standalone causal report per session. The unified
        // causal lane must emit the submission terminal before any report
        // whose pre-state includes that submission charge.
        for session in 0..MAX_SESSIONS_PER_SCOPE {
            governor
                .charge(
                    SessionId(u64::try_from(session).expect("fixture id fits")),
                    Charge::default(),
                )
                .expect("re-dirty meter");
        }
        {
            let inner = governor.inner.lock().expect("governor lock");
            assert!(matches!(
                inner.scopes["fair"].dirty_causal.first(),
                Some((_, DirtyCausalMutation::Submission(_)))
            ));
        }
        let second = governor
            .flush_scope_to_ledger(&permit, &ledger)
            .expect("causal-meter chunk");
        assert_eq!(second.appended_rows, MAX_FLUSH_ROWS);
        assert!(second.remaining_dirty);
        let third = governor
            .flush_scope_to_ledger(&permit, &ledger)
            .expect("terminal/event chunk");
        assert_eq!(third.appended_rows, 2);
        assert!(!third.remaining_dirty);
        let inner = governor.inner.lock().expect("governor lock");
        let scope = &inner.scopes["fair"];
        assert!(scope.dirty_submission_failures.is_empty());
        assert_eq!(scope.flushed_events, scope.events.len());
        assert!(scope.dirty_open_receipts.is_empty());
        assert!(scope.dirty_causal.is_empty());
    }

    #[test]
    fn dirty_open_receipt_precedes_rotated_dependent_causal_rows() {
        let governor = Governor::new();
        let permit = governor
            .open_session(test_token(20, "open-prerequisite"))
            .expect("first fixture session");
        let ledger = fs_ledger::Ledger::open(":memory:").expect("fixture ledger");
        let first = governor
            .flush_scope_to_ledger(&permit, &ledger)
            .expect("first open receipt");
        assert_eq!(first.appended_rows, 1);
        assert!(!first.remaining_dirty);
        {
            let inner = governor.inner.lock().expect("governor lock");
            assert_eq!(inner.scopes["open-prerequisite"].next_flush_lane, 1);
        }

        governor
            .open_session(test_token(21, "open-prerequisite"))
            .expect("later session with a dirty open receipt");
        for _ in 0..MAX_FLUSH_ROWS {
            governor
                .charge(SessionId(21), Charge::default())
                .expect("dependent causal report");
        }

        let second = governor
            .flush_scope_to_ledger(&permit, &ledger)
            .expect("open prerequisite plus bounded causal prefix");
        assert_eq!(second.appended_rows, MAX_FLUSH_ROWS);
        assert!(second.remaining_dirty);
        let inner = governor.inner.lock().expect("governor lock");
        let scope = &inner.scopes["open-prerequisite"];
        assert!(
            scope.dirty_open_receipts.is_empty(),
            "the later open receipt must commit in the same or an earlier transaction"
        );
        assert_eq!(
            scope.dirty_causal.len(),
            1,
            "one causal row remains because the open prerequisite occupied the first batch slot"
        );
    }

    #[test]
    fn pause_flush_defers_until_dirty_draining_submissions_commit() {
        let ledger = fs_ledger::Ledger::open(":memory:").expect("fixture ledger");
        let governor =
            super::Governor::new_durable(&ledger, DurableGovernorNonce::from_bytes([0x91; 32]))
                .expect("durable governor");
        let session = SessionId(92);
        let token = test_token(session.0, "pause-predecessor-flush");
        let open_id = governor
            .session_open_id(session, "pause-predecessor-open")
            .expect("open authority");
        let open_receipt = governor
            .open_session_gated(open_id, token, Arc::new(CancelGate::new()))
            .expect("gated session");
        let permit = open_receipt.flush_permit();
        governor
            .flush_scope_to_ledger(&permit, &ledger)
            .expect("open prerequisite");

        let submission_id = governor
            .submission_request_id(session, "draining-submission", "bounded-work")
            .expect("submission authority");
        governor
            .submit_once_durable(&ledger, submission_id, "bounded-work", Charge::default)
            .expect("completed durable submission remains dirty");
        let action_id = governor
            .pressure_action_id(session, "pause-after-submission")
            .expect("pressure authority");
        let pressure = governor
            .apply_memory_pressure(action_id, 3)
            .expect("pause request");
        let pause_request = pressure
            .events()
            .iter()
            .find_map(|event| event.pause_request_id)
            .expect("pause request authority");
        let checkpoint = solver_checkpoint(&ledger, pause_request, "checkpoint-after-bounded-work");
        governor
            .acknowledge_pause(pause_request, &ledger, &checkpoint)
            .expect("pause acknowledgement");
        let pause_authority = recovery::pause_ack_authority(pause_request);
        {
            let mut inner = governor.inner.lock().expect("governor lock");
            inner
                .scopes
                .get_mut("pause-predecessor-flush")
                .expect("fixture scope")
                .next_flush_lane = 3;
        }

        let predecessor_flush = governor
            .flush_scope_to_ledger(&permit, &ledger)
            .expect("predecessor flush");
        assert!(predecessor_flush.remaining_dirty);
        assert!(
            ledger
                .session_terminal(&submission_id.content_hash)
                .expect("verified submission terminal")
                .is_some()
        );
        assert!(
            ledger
                .session_terminal(&pause_authority)
                .expect("pause remains Pending-free and absent")
                .is_none(),
            "the control-first lane must not commit the pause ahead of its dirty predecessor"
        );

        let pause_flush = governor
            .flush_scope_to_ledger(&permit, &ledger)
            .expect("deferred pause flush");
        assert!(!pause_flush.remaining_dirty);
        assert!(
            ledger
                .session_terminal(&pause_authority)
                .expect("verified pause terminal")
                .is_some()
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)] // Directly snapshots every replay-sensitive internal counter.
    fn typed_duplicate_replay_changes_no_internal_counter_gate_or_cursor() {
        let governor = super::Governor::new();
        let session = SessionId(93);
        let open_id = governor
            .session_open_id(session, "internal-replay-open")
            .expect("open authority");
        let gate = Arc::new(CancelGate::new());
        governor
            .open_session_gated(
                open_id,
                test_token(93, "internal-replay"),
                Arc::clone(&gate),
            )
            .expect("open gated session");
        let report = governor
            .meter_report_id(session, "internal-report")
            .expect("meter authority");
        let charge = Charge {
            core_s: 2.0,
            ..Charge::default()
        };
        let meter_receipt = governor.charge(report, charge).expect("first charge");
        let before_meter_replay = {
            let inner = governor.inner.lock().expect("governor lock");
            (
                inner.next_meter_commit_ordinal,
                inner.scopes["internal-replay"].revision,
                inner.scopes["internal-replay"].dirty_causal.clone(),
                inner.meters[&session.0].snapshot(),
            )
        };
        assert_eq!(
            governor.charge(report, charge).expect("exact meter replay"),
            meter_receipt
        );
        {
            let inner = governor.inner.lock().expect("governor lock");
            assert_eq!(inner.next_meter_commit_ordinal, before_meter_replay.0);
            assert_eq!(
                inner.scopes["internal-replay"].revision,
                before_meter_replay.1
            );
            assert_eq!(
                inner.scopes["internal-replay"].dirty_causal,
                before_meter_replay.2
            );
            assert_eq!(inner.meters[&session.0].snapshot(), before_meter_replay.3);
        }

        let action = governor
            .pressure_action_id(session, "internal-pressure")
            .expect("pressure authority");
        let pressure_receipt = governor
            .apply_memory_pressure(action, 3)
            .expect("first pressure action");
        let before_pressure_replay = {
            let inner = governor.inner.lock().expect("governor lock");
            (
                inner.next_ordinal,
                inner.reserved_pause_ordinals,
                inner.scopes["internal-replay"].reserved_pause_completions,
                inner.scopes["internal-replay"].events.clone(),
                inner.scopes["internal-replay"].revision,
                Arc::clone(&inner.gates[&session.0]),
                inner.gate_generations[&session.0],
            )
        };
        assert_eq!(
            governor
                .apply_memory_pressure(action, 3)
                .expect("exact pressure replay"),
            pressure_receipt
        );
        let inner = governor.inner.lock().expect("governor lock");
        assert_eq!(inner.next_ordinal, before_pressure_replay.0);
        assert_eq!(inner.reserved_pause_ordinals, before_pressure_replay.1);
        assert_eq!(
            inner.scopes["internal-replay"].reserved_pause_completions,
            before_pressure_replay.2
        );
        assert_eq!(
            inner.scopes["internal-replay"].events,
            before_pressure_replay.3
        );
        assert_eq!(
            inner.scopes["internal-replay"].revision,
            before_pressure_replay.4
        );
        assert!(Arc::ptr_eq(
            &inner.gates[&session.0],
            &before_pressure_replay.5
        ));
        assert_eq!(inner.gate_generations[&session.0], before_pressure_replay.6);
    }

    #[test]
    #[allow(clippy::too_many_lines)] // One exact-cap matrix for both typed mutation registries.
    fn typed_report_and_action_caps_admit_replay_but_refuse_limit_plus_one() {
        let governor = super::Governor::new();
        let session = SessionId(94);
        let open_id = governor
            .session_open_id(session, "cap-open")
            .expect("open authority");
        governor
            .open_session(open_id, test_token(94, "typed-caps"))
            .expect("open session");

        let charge = Charge {
            core_s: 1.0,
            ..Charge::default()
        };
        let real_report = governor
            .meter_report_id(session, "real-report")
            .expect("real report authority");
        let meter_receipt = governor.charge(real_report, charge).expect("first report");
        let new_report = governor
            .meter_report_id(session, "limit-plus-one-report")
            .expect("uncommitted report authority");
        {
            let mut inner = governor.inner.lock().expect("governor lock");
            let open_hash = super::Governor::current_open_identity(&inner, session)
                .expect("open receipt identity");
            let reports = inner
                .meter_report_ids
                .get_mut(&session.0)
                .expect("report index");
            for ordinal in 1..MAX_METER_REPORTS_PER_SESSION {
                let hash = fs_blake3::hash_domain(
                    METER_REPORT_ID_DOMAIN,
                    &u64::try_from(ordinal)
                        .expect("bounded ordinal")
                        .to_le_bytes(),
                );
                reports.insert(MeterReportId {
                    governor_id: governor.id,
                    session,
                    session_open: open_hash,
                    generation: 0,
                    content_hash: hash,
                });
            }
            assert_eq!(reports.len(), MAX_METER_REPORTS_PER_SESSION);
        }
        assert_eq!(
            governor
                .charge(real_report, charge)
                .expect("known report replays at cap"),
            meter_receipt
        );
        assert!(matches!(
            governor.charge(new_report, charge),
            Err(SessionError::LimitExceeded {
                resource: "meter_reports_per_session",
                limit: MAX_METER_REPORTS_PER_SESSION,
                ..
            })
        ));

        let real_action = governor
            .pressure_action_id(session, "real-action")
            .expect("real action authority");
        let pressure_receipt = governor
            .apply_memory_pressure(real_action, 1)
            .expect("first pressure action");
        let new_action = governor
            .pressure_action_id(session, "limit-plus-one-action")
            .expect("uncommitted action authority");
        {
            let mut inner = governor.inner.lock().expect("governor lock");
            let open_hash = super::Governor::current_open_identity(&inner, session)
                .expect("open receipt identity");
            let actions = inner
                .pressure_action_ids
                .get_mut(&session.0)
                .expect("action index");
            for ordinal in 1..MAX_PRESSURE_ACTIONS_PER_SESSION {
                let hash = fs_blake3::hash_domain(
                    PRESSURE_ACTION_ID_DOMAIN,
                    &u64::try_from(ordinal)
                        .expect("bounded ordinal")
                        .to_le_bytes(),
                );
                actions.insert(PressureActionId {
                    governor_id: governor.id,
                    session,
                    session_open: open_hash,
                    generation: 0,
                    content_hash: hash,
                });
            }
            assert_eq!(actions.len(), MAX_PRESSURE_ACTIONS_PER_SESSION);
        }
        assert_eq!(
            governor
                .apply_memory_pressure(real_action, 1)
                .expect("known action replays at cap"),
            pressure_receipt
        );
        assert!(matches!(
            governor.apply_memory_pressure(new_action, 1),
            Err(SessionError::LimitExceeded {
                resource: "pressure_actions_per_session",
                limit: MAX_PRESSURE_ACTIONS_PER_SESSION,
                ..
            })
        ));
    }

    #[test]
    fn bounded_flush_never_splits_an_indivisible_terminal_event_group() {
        let mut terminals = Vec::new();
        let mut event_rows = MAX_FLUSH_ROWS - 3;
        let mut encoded_bytes = 0;
        assert!(
            push_bounded_terminal(
                &mut terminals,
                &mut event_rows,
                &mut encoded_bytes,
                "group-boundary",
                bounded_terminal_with_events(3),
            )
            .expect("exact event-row boundary is admitted")
        );
        assert_eq!(terminals.len(), 1);
        assert_eq!(event_rows, MAX_FLUSH_ROWS);

        let mut terminals = Vec::new();
        let mut event_rows = MAX_FLUSH_ROWS - 2;
        let mut encoded_bytes = 0;
        assert!(
            !push_bounded_terminal(
                &mut terminals,
                &mut event_rows,
                &mut encoded_bytes,
                "group-boundary",
                bounded_terminal_with_events(3),
            )
            .expect("over-bound group is deferred intact")
        );
        assert!(terminals.is_empty());
        assert_eq!(event_rows, MAX_FLUSH_ROWS - 2);
        assert_eq!(encoded_bytes, 0);
    }

    #[test]
    #[allow(clippy::too_many_lines)] // Real reopen snapshots every reservation before/after refusal.
    fn durable_pending_reopen_never_invokes_work_and_rolls_back_local_admission() {
        let path = durable_test_ledger_path("pending");
        let nonce = DurableGovernorNonce::from_bytes([0x51; 32]);
        let ledger = fs_ledger::Ledger::open(&path).expect("on-disk ledger");
        let governor = super::Governor::new_durable(&ledger, nonce).expect("durable governor");
        let session = SessionId(951);
        let token = test_token(session.0, "durable-pending");
        let open_id = governor
            .session_open_id(session, "durable-open")
            .expect("open authority");
        let open_receipt = governor
            .open_session(open_id, token.clone())
            .expect("open session");
        governor
            .flush_scope_to_ledger(&open_receipt.flush_permit(), &ledger)
            .expect("open terminal is durable before execution claims");
        let request_id = governor
            .submission_request_id(session, "pending-slot", "canonical-program")
            .expect("submission authority");
        let payload = recovery::encode_submission_payload(request_id);
        let claim = fs_ledger::session_registry::SessionMutationClaim {
            authority: request_id.content_hash,
            ledger_instance_id: ledger.checked_instance_id().expect("ledger identity"),
            governor_hash: governor.id,
            session_open_hash: request_id.session_open,
            kind: recovery::KIND_SUBMISSION,
            session: session.0,
            ledger_scope: &token.ledger_scope,
            generation: request_id.generation,
            causal_ordinal: Some(1),
            payload: &payload,
        };
        assert!(matches!(
            ledger
                .claim_session_mutation(&claim)
                .expect("fresh Pending claim"),
            fs_ledger::session_registry::SessionMutationClaimResult::Claimed { .. }
        ));
        assert_eq!(ledger.table_count("session_claims").unwrap(), 2);
        assert_eq!(ledger.table_count("session_terminals").unwrap(), 1);
        drop(governor);
        drop(ledger);

        let ledger = fs_ledger::Ledger::open(&path).expect("reopened ledger");
        let governor = super::Governor::new_durable(&ledger, nonce).expect("reopened governor");
        governor
            .recover_open(&ledger, open_id, token, None)
            .expect("recover open prerequisite");
        let before = {
            let inner = governor.inner.lock().expect("governor lock");
            (
                inner.next_submission_ordinal,
                inner.reserved_meter_ordinals,
                inner.pending_submissions[&session.0],
                inner.reserved_meter_reports[&session.0],
                inner.idempotency.len(),
                inner.idempotency_keys[&session.0].len(),
                inner.submission_admissions.len(),
                inner.retained_bytes,
                inner.scopes["durable-pending"].retained_bytes,
            )
        };
        let executions = AtomicU64::new(0);
        assert_eq!(
            governor.submit_once_durable(&ledger, request_id, "canonical-program", || {
                executions.fetch_add(1, Ordering::SeqCst);
                Charge::default()
            },),
            Err(SessionError::IndeterminateMutation {
                kind: recovery::KIND_SUBMISSION,
                authority: request_id.content_hash,
            })
        );
        assert_eq!(executions.load(Ordering::SeqCst), 0);
        let after = {
            let inner = governor.inner.lock().expect("governor lock");
            (
                inner.next_submission_ordinal,
                inner.reserved_meter_ordinals,
                inner.pending_submissions[&session.0],
                inner.reserved_meter_reports[&session.0],
                inner.idempotency.len(),
                inner.idempotency_keys[&session.0].len(),
                inner.submission_admissions.len(),
                inner.retained_bytes,
                inner.scopes["durable-pending"].retained_bytes,
            )
        };
        assert_eq!(after, before, "non-fresh claim leaves no local reservation");
        assert!(matches!(
            governor.submit_once(request_id, || {
                executions.fetch_add(1, Ordering::SeqCst);
                Charge::default()
            }),
            Err(SessionError::DurableLedgerRequired {
                kind: recovery::KIND_SUBMISSION,
                authority,
            }) if authority == request_id.content_hash
        ));
        assert_eq!(executions.load(Ordering::SeqCst), 0);
        assert_eq!(ledger.table_count("session_claims").unwrap(), 2);
        assert_eq!(ledger.table_count("session_terminals").unwrap(), 1);
        assert_eq!(ledger.table_count("events").unwrap(), 1);
    }

    #[test]
    fn durable_pending_submission_rejects_pause_flush_atomically() {
        let nonce = DurableGovernorNonce::from_bytes([0x63; 32]);
        let ledger = fs_ledger::Ledger::open(":memory:").expect("fixture ledger");
        let governor = super::Governor::new_durable(&ledger, nonce).expect("durable governor");
        let session = SessionId(953);
        let token = test_token(session.0, "durable-pending-pause");
        let open_id = governor
            .session_open_id(session, "durable-open")
            .expect("open authority");
        let open_receipt = governor
            .open_session_gated(open_id, token.clone(), Arc::new(CancelGate::new()))
            .expect("gated open");
        let permit = open_receipt.flush_permit();
        governor
            .flush_scope_to_ledger(&permit, &ledger)
            .expect("open prerequisite terminal");
        let request_id = governor
            .submission_request_id(session, "pending-slot", "pending-program")
            .expect("submission authority");
        let payload = recovery::encode_submission_payload(request_id);
        let claim = fs_ledger::session_registry::SessionMutationClaim {
            authority: request_id.content_hash,
            ledger_instance_id: ledger.checked_instance_id().expect("ledger identity"),
            governor_hash: governor.id,
            session_open_hash: request_id.session_open,
            kind: recovery::KIND_SUBMISSION,
            session: session.0,
            ledger_scope: &token.ledger_scope,
            generation: request_id.generation,
            causal_ordinal: Some(1),
            payload: &payload,
        };
        assert!(matches!(
            ledger
                .claim_session_mutation(&claim)
                .expect("crash-only Pending claim"),
            fs_ledger::session_registry::SessionMutationClaimResult::Claimed { .. }
        ));
        let action_id = governor
            .pressure_action_id(session, "pause-action")
            .expect("pressure authority");
        let pressure = governor
            .apply_memory_pressure(action_id, 3)
            .expect("pause request");
        let pause_request = pressure
            .events()
            .iter()
            .find_map(|event| event.pause_request_id)
            .expect("pause request authority");
        let checkpoint = solver_checkpoint(&ledger, pause_request, "checkpoint-after-unknown-work");
        let _acknowledgement = governor
            .acknowledge_pause(pause_request, &ledger, &checkpoint)
            .expect("fixture terminal acknowledgement");
        let before = (
            ledger.table_count("session_claims").unwrap(),
            ledger.table_count("session_terminals").unwrap(),
            ledger.table_count("session_flush_batches").unwrap(),
            ledger.table_count("events").unwrap(),
        );
        let error = governor
            .flush_scope_to_ledger(&permit, &ledger)
            .expect_err("omitted Pending work fences the complete control batch");
        assert!(error.to_string().contains("durably Pending"));
        let after = {
            let inner = governor.inner.lock().expect("governor lock");
            let scope = &inner.scopes["durable-pending-pause"];
            assert!(scope.in_flight.is_none());
            assert_eq!(scope.dirty_control.len(), 2);
            assert!(scope.dirty_causal.is_empty());
            (
                ledger.table_count("session_claims").unwrap(),
                ledger.table_count("session_terminals").unwrap(),
                ledger.table_count("session_flush_batches").unwrap(),
                ledger.table_count("events").unwrap(),
            )
        };
        assert_eq!(after, before);
    }

    #[test]
    fn future_and_truncated_terminal_codecs_fail_before_recovery_mutation() {
        let ledger = fs_ledger::Ledger::open(":memory:").expect("fixture ledger");
        let governor =
            super::Governor::new_durable(&ledger, DurableGovernorNonce::from_bytes([0x79; 32]))
                .expect("durable governor");
        let session = SessionId(954);
        let token = test_token(session.0, "codec-refusal");
        let open_id = governor
            .session_open_id(session, "codec-open")
            .expect("open authority");
        let open_receipt = governor
            .open_session(open_id, token.clone())
            .expect("open session");
        governor
            .flush_scope_to_ledger(&open_receipt.flush_permit(), &ledger)
            .expect("durable open prerequisite");
        let delta = Charge {
            core_s: 1.0,
            ..Charge::default()
        };
        let future_id = governor
            .meter_report_id(session, "future-codec")
            .expect("future authority");
        let truncated_id = governor
            .meter_report_id(session, "truncated-codec")
            .expect("truncated authority");
        for (report_id, receipt, ordinal) in [
            (future_id, 2_u32.to_le_bytes().to_vec(), 1),
            (truncated_id, vec![1, 0], 2),
        ] {
            let payload = recovery::encode_meter_payload(delta);
            let claim = fs_ledger::session_registry::SessionMutationClaim {
                authority: report_id.content_hash,
                ledger_instance_id: ledger.checked_instance_id().expect("ledger identity"),
                governor_hash: governor.id,
                session_open_hash: report_id.session_open,
                kind: recovery::KIND_METER,
                session: session.0,
                ledger_scope: &token.ledger_scope,
                generation: report_id.generation,
                causal_ordinal: Some(ordinal),
                payload: &payload,
            };
            let group = fs_ledger::session_registry::SessionTerminalGroup {
                terminal: fs_ledger::session_registry::SessionTerminalRow {
                    claim,
                    permit: None,
                    receipt: &receipt,
                },
                events: &[],
            };
            let groups = [group];
            ledger
                .append_session_terminal_batch(&fs_ledger::session_registry::SessionTerminalBatch {
                    groups: &groups,
                })
                .expect("store structurally bounded hostile codec");
        }
        let before = {
            let inner = governor.inner.lock().expect("governor lock");
            (
                inner.next_meter_commit_ordinal,
                inner.meter_reports.len(),
                inner.meter_report_ids[&session.0].len(),
                inner.meters[&session.0].snapshot(),
                inner.retained_bytes,
            )
        };
        assert_eq!(
            governor.recover_meter(&ledger, future_id, delta),
            Err(SessionError::UnsupportedTerminalSchema {
                found: 2,
                supported: recovery::TERMINAL_SCHEMA_VERSION,
            })
        );
        assert!(matches!(
            governor.recover_meter(&ledger, truncated_id, delta),
            Err(SessionError::TerminalCorrupt {
                kind: recovery::KIND_METER,
                ..
            })
        ));
        let after = {
            let inner = governor.inner.lock().expect("governor lock");
            (
                inner.next_meter_commit_ordinal,
                inner.meter_reports.len(),
                inner.meter_report_ids[&session.0].len(),
                inner.meters[&session.0].snapshot(),
                inner.retained_bytes,
            )
        };
        assert_eq!(after, before);
    }

    #[test]
    #[allow(clippy::too_many_lines)] // Real reopen covers success, failure, and commit-before-cursor replay.
    fn durable_submission_terminals_reopen_exactly_and_append_once() {
        let path = durable_test_ledger_path("terminal");
        let nonce = DurableGovernorNonce::from_bytes([0xA7; 32]);
        let ledger = fs_ledger::Ledger::open(&path).expect("on-disk ledger");
        let governor = super::Governor::new_durable(&ledger, nonce).expect("durable governor");
        let session = SessionId(952);
        let token = test_token(session.0, "durable-terminal");
        let open_id = governor
            .session_open_id(session, "durable-open")
            .expect("open authority");
        let open_receipt = governor
            .open_session(open_id, token.clone())
            .expect("open session");
        let permit = open_receipt.flush_permit();
        governor
            .flush_scope_to_ledger(&permit, &ledger)
            .expect("open prerequisite flush");

        let success_id = governor
            .submission_request_id(session, "success-slot", "success-program")
            .expect("success authority");
        let success = governor
            .submit_once_durable(&ledger, success_id, "success-program", || Charge {
                core_s: 7.0,
                mem_peak_bytes: 11,
                wall_s: 13.0,
            })
            .expect("fresh durable execution");
        let (success_receipt, meter_receipt) = match &success {
            SubmitOutcome::Executed {
                receipt,
                meter_receipt,
                ..
            } => (*receipt, meter_receipt.clone()),
            other => panic!("expected execution, got {other:?}"),
        };
        let failed_id = governor
            .submission_request_id(session, "failed-slot", "failed-program")
            .expect("failure authority");
        let failed = governor
            .submit_once_durable(&ledger, failed_id, "failed-program", || {
                panic!("durable failure evidence")
            })
            .expect("panic becomes durable terminal failure");
        let failed_receipt = match &failed {
            SubmitOutcome::Failed { receipt, .. } => *receipt,
            other => panic!("expected failure terminal, got {other:?}"),
        };
        let lane_before_flush = {
            let inner = governor.inner.lock().expect("governor lock");
            inner.scopes["durable-terminal"].next_flush_lane
        };
        let committed = governor
            .flush_scope_to_ledger(&permit, &ledger)
            .expect("submission terminal batch");
        assert_eq!(committed.appended_rows, 2);
        assert_eq!(committed.committed_terminals, 2);
        assert!(!committed.remaining_dirty);
        let counts = (
            ledger.table_count("session_claims").unwrap(),
            ledger.table_count("session_terminals").unwrap(),
            ledger.table_count("session_terminal_events").unwrap(),
            ledger.table_count("session_flush_batches").unwrap(),
            ledger.table_count("session_flush_batch_members").unwrap(),
            ledger.table_count("events").unwrap(),
        );

        // Model a process death after the atomic database commit but before
        // the in-memory cursor publication by restoring the exact prepared
        // dirty set and lane start. The identical retry must write no row.
        {
            let mut inner = governor.inner.lock().expect("governor lock");
            let scope = inner
                .scopes
                .get_mut("durable-terminal")
                .expect("fixture scope");
            scope.next_flush_lane = lane_before_flush;
            scope.dirty_causal.insert((
                meter_receipt.commit_ordinal,
                DirtyCausalMutation::Submission(success_id),
            ));
            let failed_ordinal = match failed {
                SubmitOutcome::Failed {
                    admission_ordinal, ..
                } => admission_ordinal,
                _ => unreachable!("failure checked above"),
            };
            scope
                .dirty_submission_failures
                .insert((failed_ordinal, failed_id));
        }
        let replayed_flush = governor
            .flush_scope_to_ledger(&permit, &ledger)
            .expect("commit-before-cursor retry");
        assert_eq!(replayed_flush.appended_rows, 0);
        assert_eq!(replayed_flush.committed_terminals, 2);
        assert!(!replayed_flush.remaining_dirty);
        assert_eq!(
            (
                ledger.table_count("session_claims").unwrap(),
                ledger.table_count("session_terminals").unwrap(),
                ledger.table_count("session_terminal_events").unwrap(),
                ledger.table_count("session_flush_batches").unwrap(),
                ledger.table_count("session_flush_batch_members").unwrap(),
                ledger.table_count("events").unwrap(),
            ),
            counts
        );
        drop(governor);
        drop(ledger);

        let ledger = fs_ledger::Ledger::open(&path).expect("reopened ledger");
        let governor = super::Governor::new_durable(&ledger, nonce).expect("reopened governor");
        governor
            .recover_open(&ledger, open_id, token, None)
            .expect("recover open prerequisite");
        let executions = AtomicU64::new(0);
        let success_replay = governor
            .submit_once_durable(&ledger, success_id, "success-program", || {
                executions.fetch_add(1, Ordering::SeqCst);
                Charge::default()
            })
            .expect("recover successful terminal");
        match success_replay {
            SubmitOutcome::Duplicate {
                receipt,
                meter_receipt: recovered_meter,
                ..
            } => {
                assert_eq!(receipt, success_receipt);
                assert_eq!(recovered_meter, meter_receipt);
            }
            other => panic!("expected duplicate replay, got {other:?}"),
        }
        let failed_replay = governor
            .submit_once_durable(&ledger, failed_id, "failed-program", || {
                executions.fetch_add(1, Ordering::SeqCst);
                Charge::default()
            })
            .expect("recover failed terminal");
        assert!(matches!(
            failed_replay,
            SubmitOutcome::Failed { receipt, .. } if receipt == failed_receipt
        ));
        assert_eq!(executions.load(Ordering::SeqCst), 0);

        let altered_id = governor
            .submission_request_id(session, "success-slot", "altered-program")
            .expect("same slot authority with altered payload");
        assert_eq!(altered_id.content_hash, success_id.content_hash);
        assert!(matches!(
            governor.submit_once_durable(&ledger, altered_id, "altered-program", || {
                executions.fetch_add(1, Ordering::SeqCst);
                Charge::default()
            }),
            Err(SessionError::MutationConflict {
                kind: recovery::KIND_SUBMISSION,
                ..
            })
        ));
        assert_eq!(executions.load(Ordering::SeqCst), 0);
        assert_eq!(
            (
                ledger.table_count("session_claims").unwrap(),
                ledger.table_count("session_terminals").unwrap(),
                ledger.table_count("session_terminal_events").unwrap(),
                ledger.table_count("session_flush_batches").unwrap(),
                ledger.table_count("session_flush_batch_members").unwrap(),
                ledger.table_count("events").unwrap(),
            ),
            counts
        );

        // A fresh process has no in-memory key index. Altered program bytes
        // must still conflict against the durable claim and roll back every
        // provisional local reservation before caller work can run.
        drop(governor);
        drop(ledger);
        let ledger = fs_ledger::Ledger::open(&path).expect("second reopened ledger");
        let governor = super::Governor::new_durable(&ledger, nonce).expect("second governor");
        governor
            .recover_open(
                &ledger,
                open_id,
                test_token(session.0, "durable-terminal"),
                None,
            )
            .expect("recover only the open prerequisite");
        let altered_id = governor
            .submission_request_id(session, "success-slot", "altered-program")
            .expect("same slot authority with altered durable payload");
        let before = {
            let inner = governor.inner.lock().expect("governor lock");
            (
                inner.next_submission_ordinal,
                inner.reserved_meter_ordinals,
                inner.pending_submissions[&session.0],
                inner.reserved_meter_reports[&session.0],
                inner.idempotency.len(),
                inner.idempotency_keys[&session.0].len(),
                inner.submission_admissions.len(),
                inner.retained_bytes,
                inner.scopes["durable-terminal"].retained_bytes,
            )
        };
        let executions = AtomicU64::new(0);
        assert!(matches!(
            governor.submit_once_durable(&ledger, altered_id, "altered-program", || {
                executions.fetch_add(1, Ordering::SeqCst);
                Charge::default()
            }),
            Err(SessionError::MutationConflict {
                kind: recovery::KIND_SUBMISSION,
                ..
            })
        ));
        assert_eq!(executions.load(Ordering::SeqCst), 0);
        let after = {
            let inner = governor.inner.lock().expect("governor lock");
            (
                inner.next_submission_ordinal,
                inner.reserved_meter_ordinals,
                inner.pending_submissions[&session.0],
                inner.reserved_meter_reports[&session.0],
                inner.idempotency.len(),
                inner.idempotency_keys[&session.0].len(),
                inner.submission_admissions.len(),
                inner.retained_bytes,
                inner.scopes["durable-terminal"].retained_bytes,
            )
        };
        assert_eq!(after, before);
        let failed_retained_bytes = SUBMISSION_REQUEST_RETAINED_BYTES
            + RetainedEvidence::capture("durable failure evidence")
                .preview
                .len();
        let retained_before_recovery = governor.inner.lock().expect("governor lock").retained_bytes;
        let governor = Arc::new(governor);
        let barrier = Arc::new(std::sync::Barrier::new(2));
        let mut workers = Vec::new();
        for _ in 0..2 {
            let governor = Arc::clone(&governor);
            let barrier = Arc::clone(&barrier);
            let path = path.clone();
            workers.push(std::thread::spawn(move || {
                let ledger = fs_ledger::Ledger::open(&path).expect("worker ledger handle");
                barrier.wait();
                governor.recover_submission(&ledger, failed_id, "failed-program")
            }));
        }
        for worker in workers {
            assert!(matches!(
                worker.join().expect("recovery worker joins"),
                Ok(SubmitOutcome::Failed { receipt, .. }) if receipt == failed_receipt
            ));
        }
        let inner = governor.inner.lock().expect("governor lock");
        assert_eq!(inner.idempotency.len(), 1);
        assert_eq!(
            inner.retained_bytes,
            retained_before_recovery + failed_retained_bytes,
            "concurrent recovery installs retained failure evidence once"
        );
        drop(inner);
        assert_eq!(
            (
                ledger.table_count("session_claims").unwrap(),
                ledger.table_count("session_terminals").unwrap(),
                ledger.table_count("session_terminal_events").unwrap(),
                ledger.table_count("session_flush_batches").unwrap(),
                ledger.table_count("session_flush_batch_members").unwrap(),
                ledger.table_count("events").unwrap(),
            ),
            counts
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)] // Two-session real reopen proves the governor-global fence end to end.
    fn durable_restart_fences_cross_session_work_until_global_history_is_recovered() {
        let path = durable_test_ledger_path("global-recovery-fence");
        let nonce = DurableGovernorNonce::from_bytes([0xB4; 32]);
        let ledger = fs_ledger::Ledger::open(&path).expect("on-disk ledger");
        let governor = super::Governor::new_durable(&ledger, nonce).expect("durable governor");
        let session_a = SessionId(960);
        let session_b = SessionId(961);
        let token_a = test_token(session_a.0, "global-recovery-fence");
        let token_b = test_token(session_b.0, "global-recovery-fence");
        let open_a = governor
            .session_open_id(session_a, "open-a")
            .expect("open A authority");
        let open_b = governor
            .session_open_id(session_b, "open-b")
            .expect("open B authority");
        let permit = governor
            .open_session(open_a, token_a.clone())
            .expect("open A")
            .flush_permit();
        governor
            .open_session(open_b, token_b.clone())
            .expect("open B");
        governor
            .flush_scope_to_ledger(&permit, &ledger)
            .expect("durable open prerequisites");
        let request_a = governor
            .submission_request_id(session_a, "history-a", "program-a")
            .expect("submission A authority");
        governor
            .submit_once_durable(&ledger, request_a, "program-a", || Charge {
                core_s: 1.0,
                mem_peak_bytes: 2,
                wall_s: 3.0,
            })
            .expect("historical A submission");
        governor
            .flush_scope_to_ledger(&permit, &ledger)
            .expect("historical A terminal");
        drop(governor);
        drop(ledger);

        let ledger = fs_ledger::Ledger::open(&path).expect("reopened ledger");
        let governor = super::Governor::new_durable(&ledger, nonce).expect("restarted governor");
        governor
            .recover_open(&ledger, open_b, token_b, None)
            .expect("recover only session B open");
        let request_b = governor
            .submission_request_id(session_b, "fresh-b", "program-b")
            .expect("fresh B authority");
        let executions = AtomicU64::new(0);
        assert!(matches!(
            governor.submit_once_durable(&ledger, request_b, "program-b", || {
                executions.fetch_add(1, Ordering::SeqCst);
                Charge::default()
            }),
            Err(SessionError::DurableRecoveryIncomplete {
                remaining_claims: 2
            })
        ));
        assert_eq!(executions.load(Ordering::SeqCst), 0);

        governor
            .recover_open(&ledger, open_a, token_a, None)
            .expect("recover session A open");
        governor
            .recover_submission(&ledger, request_a, "program-a")
            .expect("recover A terminal and its global meter commit");
        assert!(matches!(
            governor
                .submit_once_durable(&ledger, request_b, "program-b", || {
                    executions.fetch_add(1, Ordering::SeqCst);
                    Charge::default()
                })
                .expect("fresh B runs after complete global recovery"),
            SubmitOutcome::Executed {
                admission_ordinal: 2,
                ..
            }
        ));
        assert_eq!(executions.load(Ordering::SeqCst), 1);
    }

    #[test]
    #[allow(clippy::too_many_lines)] // Real durable setup is required before injecting the collision owner.
    fn recovered_submission_refuses_an_admission_ordinal_owned_by_another_request() {
        let path = durable_test_ledger_path("admission-collision");
        let nonce = DurableGovernorNonce::from_bytes([0xC5; 32]);
        let ledger = fs_ledger::Ledger::open(&path).expect("on-disk ledger");
        let governor = super::Governor::new_durable(&ledger, nonce).expect("durable governor");
        let session = SessionId(962);
        let token = test_token(session.0, "admission-collision");
        let open_id = governor
            .session_open_id(session, "open")
            .expect("open authority");
        let permit = governor
            .open_session(open_id, token.clone())
            .expect("open session")
            .flush_permit();
        governor
            .flush_scope_to_ledger(&permit, &ledger)
            .expect("open prerequisite");
        let historical = governor
            .submission_request_id(session, "historical", "historical-program")
            .expect("historical authority");
        assert!(matches!(
            governor
                .submit_once_durable(&ledger, historical, "historical-program", || {
                    panic!("historical failure")
                })
                .expect("panic becomes terminal"),
            SubmitOutcome::Failed {
                admission_ordinal: 1,
                ..
            }
        ));
        governor
            .flush_scope_to_ledger(&permit, &ledger)
            .expect("historical terminal");
        drop(governor);
        drop(ledger);

        let ledger = fs_ledger::Ledger::open(&path).expect("reopened ledger");
        let governor = super::Governor::new_durable(&ledger, nonce).expect("restarted governor");
        governor
            .recover_open(&ledger, open_id, token, None)
            .expect("recover open");
        let conflicting = governor
            .submission_request_id(session, "conflicting", "other-program")
            .expect("conflicting authority");
        let before = {
            let mut inner = governor.inner.lock().expect("governor lock");
            inner.submission_admissions.insert(1, conflicting);
            (
                inner.next_submission_ordinal,
                inner.idempotency.len(),
                inner.idempotency_keys[&session.0].len(),
                inner.next_meter_commit_ordinal,
                inner.retained_bytes,
                inner.recovered_durable_claims.len(),
            )
        };
        assert!(matches!(
            governor.recover_submission(&ledger, historical, "historical-program"),
            Err(SessionError::TerminalCorrupt { detail, .. })
                if detail.contains("admission ordinal 1")
        ));
        let inner = governor.inner.lock().expect("governor lock");
        assert_eq!(
            (
                inner.next_submission_ordinal,
                inner.idempotency.len(),
                inner.idempotency_keys[&session.0].len(),
                inner.next_meter_commit_ordinal,
                inner.retained_bytes,
                inner.recovered_durable_claims.len(),
            ),
            before
        );
        assert_eq!(inner.submission_admissions.get(&1), Some(&conflicting));
    }
}

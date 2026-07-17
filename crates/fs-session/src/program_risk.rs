//! Ledger-backed expansion-program risk reporting.
//!
//! The pure PR-001--PR-012 register and evaluator live in `fs-govern`.
//! This module binds one drained session snapshot to its immutable open
//! receipt, materializes both register and session report as content-addressed
//! artifacts, roots them in the lineage DAG, and publishes one typed terminal
//! plus owned alert event. The singleton report authority is content-invariant:
//! an exact retry replays, while altered report content conflicts.

use core::fmt::Write as _;

use fs_govern::program_risks::{
    AssessmentStatus, ProgramRiskAssessment, ProgramRiskId, ProgramRiskObservation,
    assess_program_risks, program_risk_register_json,
};
use fs_ledger::session_registry::{
    SessionMutationClaim, SessionTerminalBatch, SessionTerminalBatchResult, SessionTerminalGroup,
    SessionTerminalRow,
};
use fs_ledger::{EdgeRole, EventRow, ExecMode, FiveExplicits, Ledger, OpArtifactEdge, OpOutcome};

use crate::governor::{Governor, SessionOpenReceipt};
use crate::{SessionError, SessionId};

/// Artifact kind for the canonical PR-001--PR-012 register.
pub const PROGRAM_RISK_REGISTER_ARTIFACT_KIND: &str = "program-risk-register";
/// Artifact kind for one session-scoped assessment envelope.
pub const PROGRAM_RISK_SESSION_REPORT_ARTIFACT_KIND: &str = "program-risk-session-report";
/// Owned event kind for one published session report.
pub const PROGRAM_RISK_REPORT_EVENT_KIND: &str = "program-risk-report";
/// Version of the singleton report authority.
pub const PROGRAM_RISK_REPORT_IDENTITY_VERSION: u32 = 2;
/// Domain for the singleton report authority.
pub const PROGRAM_RISK_REPORT_ID_DOMAIN: &str =
    "org.frankensim.fs-session.program-risk-report-id.v2";
/// Version of the strict terminal payload/receipt codec.
pub const PROGRAM_RISK_REPORT_CODEC_VERSION: u32 = 1;
/// Version of the positional PR-001--PR-012 row order in the terminal codec.
pub const PROGRAM_RISK_REPORT_ROW_ORDER_VERSION: u32 = 1;
/// Version of the assessment-status byte tags in the terminal codec.
pub const PROGRAM_RISK_REPORT_STATUS_TAG_VERSION: u32 = 1;
/// Number of positional assessment rows carried by the v1 report codec.
pub const PROGRAM_RISK_REPORT_LOGICAL_ROWS: usize = 12;

const PROGRAM_RISK_REPORT_SLOT_TAG: u8 = 1;
const PROGRAM_RISK_REPORT_CODEC_V1_REGISTER_OFFSET: usize = 4;
const PROGRAM_RISK_REPORT_CODEC_V1_REPORT_OFFSET: usize = 36;
const PROGRAM_RISK_REPORT_CODEC_V1_LINEAGE_OFFSET: usize = 68;
const PROGRAM_RISK_REPORT_CODEC_V1_TIME_OFFSET: usize = 76;
const PROGRAM_RISK_REPORT_CODEC_V1_GENERATION_OFFSET: usize = 84;
const PROGRAM_RISK_REPORT_CODEC_V1_STATUSES_OFFSET: usize = 92;
const PROGRAM_RISK_REPORT_CODEC_V1_BYTES: usize =
    PROGRAM_RISK_REPORT_CODEC_V1_STATUSES_OFFSET + PROGRAM_RISK_REPORT_LOGICAL_ROWS;
const MAX_PROGRAM_RISK_ARTIFACT_BYTES: u64 = 1024 * 1024;
const PROGRAM_RISK_REPORT_PRODUCER_CAP: usize = 1;
const PROGRAM_RISK_LINEAGE_EDGE_CAP: usize = 2;
const PROGRAM_RISK_REGISTER_META: &str = "{\"schema\":\"frankensim.program-risk-register.v1\"}";
const PROGRAM_RISK_SESSION_REPORT_META: &str =
    "{\"schema\":\"frankensim.program-risk-session-report.v1\"}";
const PROGRAM_RISK_LINEAGE_VERSIONS_V1: &str =
    "{\"program_risk_register\":\"v1\",\"program_risk_session_report\":\"v1\"}";
const PROGRAM_RISK_LINEAGE_CAPABILITY_V1: &str = "{\"operation\":\"session.program-risk-report\"}";
const PROGRAM_RISK_V1_IDS: [ProgramRiskId; PROGRAM_RISK_REPORT_LOGICAL_ROWS] = [
    ProgramRiskId::Pr001,
    ProgramRiskId::Pr002,
    ProgramRiskId::Pr003,
    ProgramRiskId::Pr004,
    ProgramRiskId::Pr005,
    ProgramRiskId::Pr006,
    ProgramRiskId::Pr007,
    ProgramRiskId::Pr008,
    ProgramRiskId::Pr009,
    ProgramRiskId::Pr010,
    ProgramRiskId::Pr011,
    ProgramRiskId::Pr012,
];

#[derive(Debug, Clone, Copy)]
struct ProgramRiskReportIdIdentitySource {
    governor_id: fs_blake3::ContentHash,
    session: SessionId,
    session_open: fs_blake3::ContentHash,
}

fn program_risk_report_authority_identity(
    source: &ProgramRiskReportIdIdentitySource,
) -> fs_blake3::ContentHash {
    let mut canonical = Vec::with_capacity(32 + 8 + 32 + 1);
    canonical.extend_from_slice(source.governor_id.as_bytes());
    canonical.extend_from_slice(&source.session.0.to_le_bytes());
    canonical.extend_from_slice(source.session_open.as_bytes());
    canonical.push(PROGRAM_RISK_REPORT_SLOT_TAG);
    fs_blake3::hash_domain(PROGRAM_RISK_REPORT_ID_DOMAIN, &canonical)
}

#[allow(dead_code)]
fn classify_program_risk_report_id_identity_fields(source: &ProgramRiskReportIdIdentitySource) {
    let ProgramRiskReportIdIdentitySource {
        governor_id,
        session,
        session_open,
    } = source;
    let _ = (governor_id, session, session_open);
}

fn program_risk_report_identity_transport_is_current(id: ProgramRiskReportId) -> bool {
    id.content_hash
        == program_risk_report_authority_identity(&ProgramRiskReportIdIdentitySource {
            governor_id: id.governor_id,
            session: id.session,
            session_open: id.session_open,
        })
}

/// Owner-local identity declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const PROGRAM_RISK_REPORT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-session:program-risk-report-id",
    "version_const=PROGRAM_RISK_REPORT_IDENTITY_VERSION",
    "version=2",
    "domain=org.frankensim.fs-session.program-risk-report-id.v2",
    "domain_const=PROGRAM_RISK_REPORT_ID_DOMAIN",
    "encoder=program_risk_report_authority_identity",
    "encoder_helpers=none",
    "schema_functions=Governor::write_program_risk_session_end_report,Governor::recover_program_risk_report,crates/fs-blake3/src/lib.rs#hash_domain",
    "schema_constants=PROGRAM_RISK_REPORT_IDENTITY_VERSION,PROGRAM_RISK_REPORT_ID_DOMAIN,PROGRAM_RISK_REPORT_SLOT_TAG,PROGRAM_RISK_REPORT_CODEC_VERSION,PROGRAM_RISK_REPORT_ROW_ORDER_VERSION,PROGRAM_RISK_REPORT_STATUS_TAG_VERSION,PROGRAM_RISK_REPORT_LOGICAL_ROWS",
    "schema_dependencies=fs-session:durable-governor-id,fs-session:session-open-receipt",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=ProgramRiskReportIdIdentitySource",
    "source_fields=ProgramRiskReportIdIdentitySource.governor_id:semantic,ProgramRiskReportIdIdentitySource.session:semantic,ProgramRiskReportIdIdentitySource.session_open:semantic",
    "source_bindings=ProgramRiskReportIdIdentitySource.governor_id>governor-id,ProgramRiskReportIdIdentitySource.session>session,ProgramRiskReportIdIdentitySource.session_open>session-open",
    "external_semantic_fields=identity-version,digest-domain,canonical-field-order,singleton-slot",
    "semantic_fields=identity-version,digest-domain,canonical-field-order,singleton-slot,governor-id,session,session-open",
    "excluded_fields=none",
    "consumers=Governor::write_program_risk_session_end_report,Governor::recover_program_risk_report,ProgramRiskReportReceipt,fs-ledger::session_registry::SessionMutationClaim",
    "mutations=identity-version:crates/fs-session/src/program_risk.rs#program_risk_report_identity_fields_move_independently,digest-domain:crates/fs-session/src/program_risk.rs#program_risk_report_identity_fields_move_independently,canonical-field-order:crates/fs-session/src/program_risk.rs#program_risk_report_identity_fields_move_independently,singleton-slot:crates/fs-session/src/program_risk.rs#program_risk_report_identity_fields_move_independently,governor-id:crates/fs-session/src/program_risk.rs#program_risk_report_identity_fields_move_independently,session:crates/fs-session/src/program_risk.rs#program_risk_report_identity_fields_move_independently,session-open:crates/fs-session/src/program_risk.rs#program_risk_report_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_program_risk_report_id_identity_fields",
    "transport_guard=program_risk_report_identity_transport_is_current",
    "version_guard=crates/fs-session/src/program_risk.rs#program_risk_report_identity_version_and_transport_fail_closed",
    "coupling_surface=fs-session:program-risk-report-id",
];

/// Opaque singleton authority for one session's program-risk snapshot slot.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ProgramRiskReportId {
    governor_id: fs_blake3::ContentHash,
    session: SessionId,
    session_open: fs_blake3::ContentHash,
    content_hash: fs_blake3::ContentHash,
}

impl ProgramRiskReportId {
    /// Session whose report slot this authority names.
    #[must_use]
    pub const fn session(self) -> SessionId {
        self.session
    }

    /// Domain-separated singleton authority.
    #[must_use]
    pub const fn content_hash(self) -> fs_blake3::ContentHash {
        self.content_hash
    }
}

/// One automatically surfaced non-green program-risk row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ProgramRiskAlert {
    /// Canonical risk id.
    pub id: ProgramRiskId,
    /// Exact fail-closed disposition.
    pub status: AssessmentStatus,
}

/// Semantic receipt for one persisted session report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProgramRiskReportReceipt {
    report_id: ProgramRiskReportId,
    register_artifact: fs_blake3::ContentHash,
    report_artifact: fs_blake3::ContentHash,
    lineage_op: i64,
    logical_time: i64,
    generation: u64,
    statuses: [AssessmentStatus; PROGRAM_RISK_REPORT_LOGICAL_ROWS],
}

impl ProgramRiskReportReceipt {
    /// Singleton report authority.
    #[must_use]
    pub const fn report_id(&self) -> ProgramRiskReportId {
        self.report_id
    }

    /// Content address of the canonical twelve-row register.
    #[must_use]
    pub const fn register_artifact(&self) -> fs_blake3::ContentHash {
        self.register_artifact
    }

    /// Content address of the session-scoped assessment envelope.
    #[must_use]
    pub const fn report_artifact(&self) -> fs_blake3::ContentHash {
        self.report_artifact
    }

    /// Ledger operation that roots register-input and report-output lineage.
    #[must_use]
    pub const fn lineage_op(&self) -> i64 {
        self.lineage_op
    }

    /// Explicit deterministic event time supplied by the caller.
    #[must_use]
    pub const fn logical_time(&self) -> i64 {
        self.logical_time
    }

    /// Execution generation in which the immutable snapshot was published.
    #[must_use]
    pub const fn generation(&self) -> u64 {
        self.generation
    }

    /// Every non-green row, in canonical risk order.
    #[must_use]
    pub fn alerts(&self) -> Vec<ProgramRiskAlert> {
        PROGRAM_RISK_V1_IDS
            .into_iter()
            .zip(self.statuses)
            .filter_map(|(id, status)| status.is_alert().then_some(ProgramRiskAlert { id, status }))
            .collect()
    }

    /// Number of surfaced non-green rows.
    #[must_use]
    pub fn alert_count(&self) -> usize {
        self.statuses
            .iter()
            .filter(|status| status.is_alert())
            .count()
    }
}

/// Storage result for a report publication attempt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProgramRiskReportDisposition {
    /// The terminal, event, and batch marker were newly committed.
    Committed,
    /// The exact terminal batch already existed; zero rows were appended.
    Replayed,
}

/// Receipt plus the storage observation for one write.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProgramRiskReportWrite {
    /// Semantic report receipt.
    pub receipt: ProgramRiskReportReceipt,
    /// Whether this call committed or replayed the immutable batch.
    pub disposition: ProgramRiskReportDisposition,
}

fn status_tag(status: AssessmentStatus) -> u8 {
    match status {
        AssessmentStatus::Clear => 0,
        AssessmentStatus::Triggered => 1,
        AssessmentStatus::Missing => 2,
        AssessmentStatus::Duplicate => 3,
        AssessmentStatus::NonFinite => 4,
        AssessmentStatus::UnitMismatch => 5,
        AssessmentStatus::OutOfRange => 6,
        AssessmentStatus::UnderSampled => 7,
    }
}

fn status_code_v1(status: AssessmentStatus) -> &'static str {
    match status {
        AssessmentStatus::Clear => "clear",
        AssessmentStatus::Triggered => "triggered",
        AssessmentStatus::Missing => "missing",
        AssessmentStatus::Duplicate => "duplicate",
        AssessmentStatus::NonFinite => "non-finite",
        AssessmentStatus::UnitMismatch => "unit-mismatch",
        AssessmentStatus::OutOfRange => "out-of-range",
        AssessmentStatus::UnderSampled => "under-sampled",
    }
}

fn status_from_tag(
    tag: u8,
    authority: fs_blake3::ContentHash,
) -> Result<AssessmentStatus, SessionError> {
    match tag {
        0 => Ok(AssessmentStatus::Clear),
        1 => Ok(AssessmentStatus::Triggered),
        2 => Ok(AssessmentStatus::Missing),
        3 => Ok(AssessmentStatus::Duplicate),
        4 => Ok(AssessmentStatus::NonFinite),
        5 => Ok(AssessmentStatus::UnitMismatch),
        6 => Ok(AssessmentStatus::OutOfRange),
        7 => Ok(AssessmentStatus::UnderSampled),
        other => Err(SessionError::TerminalCorrupt {
            kind: crate::governor::recovery::KIND_PROGRAM_RISK_REPORT,
            authority,
            detail: format!("unknown program-risk assessment status tag {other}"),
        }),
    }
}

fn assessment_statuses(
    assessment: &ProgramRiskAssessment,
) -> Result<[AssessmentStatus; PROGRAM_RISK_REPORT_LOGICAL_ROWS], SessionError> {
    assessment
        .rows()
        .iter()
        .map(|row| row.status)
        .collect::<Vec<_>>()
        .try_into()
        .map_err(|rows: Vec<AssessmentStatus>| SessionError::Persistence {
            what: format!(
                "program-risk evaluator returned {} rows; expected {PROGRAM_RISK_REPORT_LOGICAL_ROWS}",
                rows.len()
            ),
        })
}

fn encode_terminal_receipt(receipt: &ProgramRiskReportReceipt) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(PROGRAM_RISK_REPORT_CODEC_V1_BYTES);
    bytes.extend_from_slice(&PROGRAM_RISK_REPORT_CODEC_VERSION.to_le_bytes());
    bytes.extend_from_slice(receipt.register_artifact.as_bytes());
    bytes.extend_from_slice(receipt.report_artifact.as_bytes());
    bytes.extend_from_slice(&receipt.lineage_op.to_le_bytes());
    bytes.extend_from_slice(&receipt.logical_time.to_le_bytes());
    bytes.extend_from_slice(&receipt.generation.to_le_bytes());
    bytes.extend(receipt.statuses.into_iter().map(status_tag));
    bytes
}

fn decode_terminal_receipt(
    report_id: ProgramRiskReportId,
    bytes: &[u8],
) -> Result<ProgramRiskReportReceipt, SessionError> {
    let corrupt = |detail: String| SessionError::TerminalCorrupt {
        kind: crate::governor::recovery::KIND_PROGRAM_RISK_REPORT,
        authority: report_id.content_hash,
        detail,
    };
    if bytes.len() < core::mem::size_of::<u32>() {
        return Err(corrupt(format!(
            "program-risk terminal has {} bytes; a codec version needs at least {}",
            bytes.len(),
            core::mem::size_of::<u32>()
        )));
    }
    let version = u32::from_le_bytes(bytes[0..4].try_into().expect("checked codec length"));
    match version {
        1 => decode_terminal_receipt_v1(report_id, bytes),
        found if found > PROGRAM_RISK_REPORT_CODEC_VERSION => {
            Err(SessionError::UnsupportedTerminalSchema {
                found,
                supported: PROGRAM_RISK_REPORT_CODEC_VERSION,
            })
        }
        found => Err(corrupt(format!(
            "program-risk terminal codec v{found} is not supported"
        ))),
    }
}

fn decode_terminal_receipt_v1(
    report_id: ProgramRiskReportId,
    bytes: &[u8],
) -> Result<ProgramRiskReportReceipt, SessionError> {
    let corrupt = |detail: String| SessionError::TerminalCorrupt {
        kind: crate::governor::recovery::KIND_PROGRAM_RISK_REPORT,
        authority: report_id.content_hash,
        detail,
    };
    if bytes.len() != PROGRAM_RISK_REPORT_CODEC_V1_BYTES {
        return Err(corrupt(format!(
            "program-risk terminal codec v1 has {} bytes; expected {PROGRAM_RISK_REPORT_CODEC_V1_BYTES}",
            bytes.len()
        )));
    }
    let version = u32::from_le_bytes(bytes[0..4].try_into().expect("checked codec length"));
    if version != 1 {
        return Err(corrupt(format!(
            "program-risk terminal v1 decoder received codec v{version}"
        )));
    }
    let register_artifact = fs_blake3::ContentHash(
        bytes[PROGRAM_RISK_REPORT_CODEC_V1_REGISTER_OFFSET
            ..PROGRAM_RISK_REPORT_CODEC_V1_REPORT_OFFSET]
            .try_into()
            .expect("checked register-hash width"),
    );
    let report_artifact = fs_blake3::ContentHash(
        bytes[PROGRAM_RISK_REPORT_CODEC_V1_REPORT_OFFSET
            ..PROGRAM_RISK_REPORT_CODEC_V1_LINEAGE_OFFSET]
            .try_into()
            .expect("checked report-hash width"),
    );
    let lineage_op = i64::from_le_bytes(
        bytes
            [PROGRAM_RISK_REPORT_CODEC_V1_LINEAGE_OFFSET..PROGRAM_RISK_REPORT_CODEC_V1_TIME_OFFSET]
            .try_into()
            .expect("checked lineage-op width"),
    );
    if lineage_op <= 0 {
        return Err(corrupt(format!(
            "program-risk lineage op must be positive, found {lineage_op}"
        )));
    }
    let logical_time = i64::from_le_bytes(
        bytes[PROGRAM_RISK_REPORT_CODEC_V1_TIME_OFFSET
            ..PROGRAM_RISK_REPORT_CODEC_V1_GENERATION_OFFSET]
            .try_into()
            .expect("checked logical-time width"),
    );
    let generation = u64::from_le_bytes(
        bytes[PROGRAM_RISK_REPORT_CODEC_V1_GENERATION_OFFSET
            ..PROGRAM_RISK_REPORT_CODEC_V1_STATUSES_OFFSET]
            .try_into()
            .expect("checked generation width"),
    );
    let mut statuses = [AssessmentStatus::Missing; PROGRAM_RISK_REPORT_LOGICAL_ROWS];
    for (slot, tag) in statuses.iter_mut().zip(
        &bytes[PROGRAM_RISK_REPORT_CODEC_V1_STATUSES_OFFSET..PROGRAM_RISK_REPORT_CODEC_V1_BYTES],
    ) {
        *slot = status_from_tag(*tag, report_id.content_hash)?;
    }
    Ok(ProgramRiskReportReceipt {
        report_id,
        register_artifact,
        report_artifact,
        lineage_op,
        logical_time,
        generation,
        statuses,
    })
}

fn report_id(governor: &Governor, open_receipt: &SessionOpenReceipt) -> ProgramRiskReportId {
    let source = ProgramRiskReportIdIdentitySource {
        governor_id: governor.identity(),
        session: open_receipt.open_id().session(),
        session_open: open_receipt.content_hash(),
    };
    ProgramRiskReportId {
        governor_id: source.governor_id,
        session: source.session,
        session_open: source.session_open,
        content_hash: program_risk_report_authority_identity(&source),
    }
}

fn report_artifact_json(
    session: SessionId,
    ledger_scope: &str,
    open_receipt: &SessionOpenReceipt,
    register_artifact: fs_blake3::ContentHash,
    logical_time: i64,
    generation: u64,
    assessment: &ProgramRiskAssessment,
) -> String {
    format!(
        "{{\"schema\":\"frankensim.program-risk-session-report.v1\",\"session\":{},\"ledger_scope\":\"{}\",\"session_open\":\"{}\",\"register_artifact\":\"{}\",\"logical_time\":{},\"generation\":{},\"assessment\":{}}}",
        session.0,
        json_escape(ledger_scope),
        open_receipt.content_hash(),
        register_artifact,
        logical_time,
        generation,
        assessment.to_json(),
    )
}

fn json_escape(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{0000}'..='\u{001f}' => {
                write!(out, "\\u{:04x}", u32::from(character))
                    .expect("writing to a String is infallible");
            }
            _ => out.push(character),
        }
    }
    out
}

fn lineage_ir(
    report_id: ProgramRiskReportId,
    register_artifact: fs_blake3::ContentHash,
    report_artifact: fs_blake3::ContentHash,
) -> String {
    format!(
        "{{\"op\":\"session.program-risk-report\",\"authority\":\"{}\",\"session\":{},\"register_artifact\":\"{}\",\"report_artifact\":\"{}\"}}",
        report_id.content_hash, report_id.session.0, register_artifact, report_artifact,
    )
}

fn lineage_budget(register_len: usize, report_len: usize) -> String {
    format!("{{\"register_bytes\":{register_len},\"report_bytes\":{report_len}}}")
}

fn map_ledger(context: &str, error: fs_ledger::LedgerError) -> SessionError {
    SessionError::Persistence {
        what: format!("{context}: {error}"),
    }
}

fn verify_artifact(
    ledger: &Ledger,
    hash: fs_blake3::ContentHash,
    kind: &str,
    expected: &[u8],
) -> Result<(), SessionError> {
    let expected_len = u64::try_from(expected.len()).map_err(|_| SessionError::Persistence {
        what: "expected program-risk artifact length does not fit u64".to_string(),
    })?;
    let info = ledger
        .artifact_info(&hash)
        .map_err(|error| map_ledger("program-risk artifact metadata", error))?
        .ok_or_else(|| SessionError::TerminalCorrupt {
            kind: crate::governor::recovery::KIND_PROGRAM_RISK_REPORT,
            authority: hash,
            detail: "referenced program-risk artifact is missing".to_string(),
        })?;
    if info.kind != kind || info.len != expected_len {
        return Err(SessionError::TerminalCorrupt {
            kind: crate::governor::recovery::KIND_PROGRAM_RISK_REPORT,
            authority: hash,
            detail: format!(
                "artifact envelope is kind {:?} / {} bytes; expected {kind:?} / {} bytes",
                info.kind,
                info.len,
                expected.len()
            ),
        });
    }
    let stored = ledger
        .get_artifact_bounded(&hash, expected_len)
        .map_err(|error| map_ledger("program-risk artifact read", error))?
        .ok_or_else(|| SessionError::TerminalCorrupt {
            kind: crate::governor::recovery::KIND_PROGRAM_RISK_REPORT,
            authority: hash,
            detail: "referenced program-risk artifact disappeared".to_string(),
        })?;
    if stored != expected {
        return Err(SessionError::TerminalCorrupt {
            kind: crate::governor::recovery::KIND_PROGRAM_RISK_REPORT,
            authority: hash,
            detail: "program-risk artifact bytes do not match their expected canonical form"
                .to_string(),
        });
    }
    Ok(())
}

fn read_historic_artifact(
    ledger: &Ledger,
    hash: fs_blake3::ContentHash,
    expected_kind: &str,
) -> Result<Vec<u8>, SessionError> {
    let info = ledger
        .artifact_info(&hash)
        .map_err(|error| map_ledger("historic program-risk artifact metadata", error))?
        .ok_or_else(|| SessionError::TerminalCorrupt {
            kind: crate::governor::recovery::KIND_PROGRAM_RISK_REPORT,
            authority: hash,
            detail: "referenced historic program-risk artifact is missing".to_string(),
        })?;
    if info.kind != expected_kind {
        return Err(SessionError::TerminalCorrupt {
            kind: crate::governor::recovery::KIND_PROGRAM_RISK_REPORT,
            authority: hash,
            detail: format!(
                "historic artifact kind is {:?}; expected {expected_kind:?}",
                info.kind
            ),
        });
    }
    if info.len > MAX_PROGRAM_RISK_ARTIFACT_BYTES {
        return Err(SessionError::TerminalCorrupt {
            kind: crate::governor::recovery::KIND_PROGRAM_RISK_REPORT,
            authority: hash,
            detail: format!(
                "historic artifact has {} bytes; limit is {MAX_PROGRAM_RISK_ARTIFACT_BYTES}",
                info.len
            ),
        });
    }
    let artifact_len = usize::try_from(info.len).map_err(|_| SessionError::TerminalCorrupt {
        kind: crate::governor::recovery::KIND_PROGRAM_RISK_REPORT,
        authority: hash,
        detail: format!("historic artifact length {} does not fit usize", info.len),
    })?;
    let stored = ledger
        .get_artifact_bounded(&hash, MAX_PROGRAM_RISK_ARTIFACT_BYTES)
        .map_err(|error| map_ledger("historic program-risk artifact read", error))?
        .ok_or_else(|| SessionError::TerminalCorrupt {
            kind: crate::governor::recovery::KIND_PROGRAM_RISK_REPORT,
            authority: hash,
            detail: "referenced historic program-risk artifact disappeared".to_string(),
        })?;
    if stored.len() != artifact_len || fs_blake3::hash_bytes(&stored) != hash {
        return Err(SessionError::TerminalCorrupt {
            kind: crate::governor::recovery::KIND_PROGRAM_RISK_REPORT,
            authority: hash,
            detail: "historic artifact length or content address disagrees with its envelope"
                .to_string(),
        });
    }
    Ok(stored)
}

const PROGRAM_RISK_V1_ROW_CODES: [&str; PROGRAM_RISK_REPORT_LOGICAL_ROWS] = [
    "PR-001", "PR-002", "PR-003", "PR-004", "PR-005", "PR-006", "PR-007", "PR-008", "PR-009",
    "PR-010", "PR-011", "PR-012",
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum HistoricComparisonV1 {
    GreaterThanOrEqual,
    LessThan,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum HistoricDomainV1 {
    NonNegativeReal,
    Fraction,
    NonNegativeInteger,
}

impl HistoricDomainV1 {
    fn admits(self, value: f64) -> bool {
        match self {
            Self::NonNegativeReal => value >= 0.0,
            Self::Fraction => (0.0..=1.0).contains(&value),
            Self::NonNegativeInteger => value >= 0.0 && value.fract() == 0.0,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct HistoricTriggerV1 {
    comparison: HistoricComparisonV1,
    threshold_bits: u64,
    unit: String,
    domain: HistoricDomainV1,
    min_samples: u64,
}

impl HistoricTriggerV1 {
    fn threshold(&self) -> f64 {
        f64::from_bits(self.threshold_bits)
    }

    fn assess(&self, value: f64, unit: &str, samples: u64) -> AssessmentStatus {
        if unit != self.unit {
            AssessmentStatus::UnitMismatch
        } else if !self.domain.admits(value) {
            AssessmentStatus::OutOfRange
        } else if samples < self.min_samples {
            AssessmentStatus::UnderSampled
        } else {
            let triggered = match self.comparison {
                HistoricComparisonV1::GreaterThanOrEqual => value >= self.threshold(),
                HistoricComparisonV1::LessThan => value < self.threshold(),
            };
            if triggered {
                AssessmentStatus::Triggered
            } else {
                AssessmentStatus::Clear
            }
        }
    }
}

struct HistoricRegisterV1 {
    triggers: Vec<HistoricTriggerV1>,
}

struct HistoricReportExpected<'a> {
    session: SessionId,
    ledger_scope: &'a str,
    session_open: fs_blake3::ContentHash,
    register_artifact: fs_blake3::ContentHash,
    logical_time: i64,
    generation: u64,
    statuses: &'a [AssessmentStatus; PROGRAM_RISK_REPORT_LOGICAL_ROWS],
}

/// Cursor for the exact whitespace-free JSON emitted by the v1 writers.
/// Historic decoders stay frozen; a future schema gets a separate decoder.
struct ProgramRiskJsonCursor<'a> {
    input: &'a str,
    offset: usize,
}

impl ProgramRiskJsonCursor<'_> {
    fn take(&mut self, expected: &str) -> Option<()> {
        self.input
            .get(self.offset..)?
            .starts_with(expected)
            .then(|| self.offset += expected.len())
    }

    fn is_finished(&self) -> bool {
        self.offset == self.input.len()
    }

    fn canonical_string(&mut self) -> Option<String> {
        let start = self.offset;
        self.take("\"")?;
        let mut value = String::new();
        loop {
            let rest = self.input.get(self.offset..)?;
            let character = rest.chars().next()?;
            match character {
                '"' => {
                    self.offset += 1;
                    break;
                }
                '\\' => {
                    self.offset += 1;
                    let escape = self.input.get(self.offset..)?.chars().next()?;
                    self.offset += escape.len_utf8();
                    match escape {
                        '"' => value.push('"'),
                        '\\' => value.push('\\'),
                        'n' => value.push('\n'),
                        'r' => value.push('\r'),
                        't' => value.push('\t'),
                        'u' => {
                            let scalar = u32::from(self.take_hex_quad()?);
                            if scalar > 0x1f {
                                return None;
                            }
                            value.push(char::from_u32(scalar)?);
                        }
                        _ => return None,
                    }
                }
                '\u{0000}'..='\u{001f}' => return None,
                _ => {
                    value.push(character);
                    self.offset += character.len_utf8();
                }
            }
        }
        let consumed = self.input.get(start..self.offset)?;
        let canonical = format!("\"{}\"", json_escape(&value));
        (consumed == canonical).then_some(value)
    }

    fn take_hex_quad(&mut self) -> Option<u16> {
        let end = self.offset.checked_add(4)?;
        let hex = self.input.get(self.offset..end)?;
        if !hex
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
        {
            return None;
        }
        self.offset = end;
        u16::from_str_radix(hex, 16).ok()
    }

    fn canonical_u64(&mut self) -> Option<u64> {
        let start = self.offset;
        while self
            .input
            .as_bytes()
            .get(self.offset)
            .is_some_and(u8::is_ascii_digit)
        {
            self.offset += 1;
        }
        let digits = self.input.get(start..self.offset)?;
        let value = digits.parse::<u64>().ok()?;
        (value.to_string() == digits).then_some(value)
    }

    fn canonical_i64(&mut self) -> Option<i64> {
        let start = self.offset;
        if self.input.as_bytes().get(self.offset) == Some(&b'-') {
            self.offset += 1;
        }
        while self
            .input
            .as_bytes()
            .get(self.offset)
            .is_some_and(u8::is_ascii_digit)
        {
            self.offset += 1;
        }
        let digits = self.input.get(start..self.offset)?;
        let value = digits.parse::<i64>().ok()?;
        (value.to_string() == digits).then_some(value)
    }

    fn canonical_f64(&mut self) -> Option<f64> {
        let start = self.offset;
        if self.input.as_bytes().get(self.offset) == Some(&b'-') {
            self.offset += 1;
        }
        match self.input.as_bytes().get(self.offset)? {
            b'0' => self.offset += 1,
            b'1'..=b'9' => {
                self.offset += 1;
                while self
                    .input
                    .as_bytes()
                    .get(self.offset)
                    .is_some_and(u8::is_ascii_digit)
                {
                    self.offset += 1;
                }
            }
            _ => return None,
        }
        if self.input.as_bytes().get(self.offset) == Some(&b'.') {
            self.offset += 1;
            let fraction_start = self.offset;
            while self
                .input
                .as_bytes()
                .get(self.offset)
                .is_some_and(u8::is_ascii_digit)
            {
                self.offset += 1;
            }
            if self.offset == fraction_start {
                return None;
            }
        }
        if matches!(self.input.as_bytes().get(self.offset), Some(b'e' | b'E')) {
            self.offset += 1;
            if matches!(self.input.as_bytes().get(self.offset), Some(b'+' | b'-')) {
                self.offset += 1;
            }
            let exponent_start = self.offset;
            while self
                .input
                .as_bytes()
                .get(self.offset)
                .is_some_and(u8::is_ascii_digit)
            {
                self.offset += 1;
            }
            if self.offset == exponent_start {
                return None;
            }
        }
        let token = self.input.get(start..self.offset)?;
        let value = token.parse::<f64>().ok()?;
        (value.is_finite() && value.to_string() == token).then_some(value)
    }

    fn canonical_optional_f64(&mut self) -> Option<Option<f64>> {
        if self.input.get(self.offset..)?.starts_with("null") {
            self.take("null")?;
            Some(None)
        } else {
            self.canonical_f64().map(Some)
        }
    }

    fn canonical_optional_u64(&mut self) -> Option<Option<u64>> {
        if self.input.get(self.offset..)?.starts_with("null") {
            self.take("null")?;
            Some(None)
        } else {
            self.canonical_u64().map(Some)
        }
    }

    fn canonical_optional_string(&mut self) -> Option<Option<String>> {
        if self.input.get(self.offset..)?.starts_with("null") {
            self.take("null")?;
            Some(None)
        } else {
            self.canonical_string().map(Some)
        }
    }

    fn canonical_bool(&mut self) -> Option<bool> {
        if self.input.get(self.offset..)?.starts_with("true") {
            self.take("true")?;
            Some(true)
        } else {
            self.take("false")?;
            Some(false)
        }
    }
}

fn historic_corrupt(authority: fs_blake3::ContentHash, detail: impl Into<String>) -> SessionError {
    SessionError::TerminalCorrupt {
        kind: crate::governor::recovery::KIND_PROGRAM_RISK_REPORT,
        authority,
        detail: detail.into(),
    }
}

fn bounded_utf8_prefix(value: &str, max_bytes: usize) -> &str {
    let mut end = 0;
    for (offset, character) in value.char_indices() {
        let next = offset + character.len_utf8();
        if next > max_bytes {
            break;
        }
        end = next;
    }
    &value[..end]
}

fn parse_historic_trigger_v1(parser: &mut ProgramRiskJsonCursor<'_>) -> Option<HistoricTriggerV1> {
    parser.take("{\"comparison\":")?;
    let comparison = match parser.canonical_string()?.as_str() {
        ">=" => HistoricComparisonV1::GreaterThanOrEqual,
        "<" => HistoricComparisonV1::LessThan,
        _ => return None,
    };
    parser.take(",\"threshold\":")?;
    let threshold = parser.canonical_f64()?;
    parser.take(",\"unit\":")?;
    let unit = parser.canonical_string()?;
    parser.take(",\"domain\":")?;
    let domain = match parser.canonical_string()?.as_str() {
        "non-negative-real" => HistoricDomainV1::NonNegativeReal,
        "fraction-0-to-1" => HistoricDomainV1::Fraction,
        "non-negative-integer" => HistoricDomainV1::NonNegativeInteger,
        _ => return None,
    };
    parser.take(",\"min_samples\":")?;
    let min_samples = parser.canonical_u64()?;
    parser.take("}")?;
    if unit.trim().is_empty()
        || unit.len() > fs_govern::program_risks::MAX_OBSERVED_UNIT_PREVIEW_BYTES
        || min_samples == 0
        || !domain.admits(threshold)
    {
        return None;
    }
    Some(HistoricTriggerV1 {
        comparison,
        threshold_bits: threshold.to_bits(),
        unit,
        domain,
        min_samples,
    })
}

fn decode_historic_register(
    bytes: &[u8],
    authority: fs_blake3::ContentHash,
) -> Result<HistoricRegisterV1, SessionError> {
    let input = core::str::from_utf8(bytes)
        .map_err(|_| historic_corrupt(authority, "historic program-risk register is not UTF-8"))?;
    let mut parser = ProgramRiskJsonCursor { input, offset: 0 };
    parser
        .take("{\"schema\":")
        .ok_or_else(|| historic_corrupt(authority, "historic register lacks a canonical schema"))?;
    let schema = parser.canonical_string().ok_or_else(|| {
        historic_corrupt(
            authority,
            "historic register schema is not a canonical string",
        )
    })?;
    if schema != "frankensim.program-risk-register.v1" {
        return Err(historic_corrupt(
            authority,
            format!(
                "unsupported historic program-risk register schema {:?} ({} bytes)",
                bounded_utf8_prefix(&schema, 64),
                schema.len()
            ),
        ));
    }
    decode_historic_register_v1(&mut parser).ok_or_else(|| {
        historic_corrupt(
            authority,
            "historic program-risk register v1 violates its frozen canonical schema",
        )
    })
}

fn decode_historic_register_v1(
    parser: &mut ProgramRiskJsonCursor<'_>,
) -> Option<HistoricRegisterV1> {
    parser.take(",\"risks\":[")?;
    let mut triggers = Vec::with_capacity(PROGRAM_RISK_REPORT_LOGICAL_ROWS);
    for (index, expected_id) in PROGRAM_RISK_V1_ROW_CODES.into_iter().enumerate() {
        if index > 0 {
            parser.take(",")?;
        }
        parser.take("{\"id\":")?;
        if parser.canonical_string()? != expected_id {
            return None;
        }
        parser.take(",\"name\":")?;
        if parser.canonical_string()?.trim().is_empty() {
            return None;
        }
        parser.take(",\"owner\":{\"role\":")?;
        if parser.canonical_string()?.trim().is_empty() {
            return None;
        }
        parser.take(",\"bead_id\":")?;
        if !parser.canonical_string()?.starts_with("frankensim-") {
            return None;
        }
        parser.take("},\"likelihood\":")?;
        if !(1..=5).contains(&parser.canonical_u64()?) {
            return None;
        }
        parser.take(",\"impact\":")?;
        if !(1..=5).contains(&parser.canonical_u64()?) {
            return None;
        }
        parser.take(",\"leading_indicator\":")?;
        if parser.canonical_string()?.trim().is_empty() {
            return None;
        }
        parser.take(",\"trigger\":")?;
        triggers.push(parse_historic_trigger_v1(parser)?);
        parser.take(",\"mitigation\":")?;
        if parser.canonical_string()?.trim().is_empty() {
            return None;
        }
        parser.take(",\"contingency\":")?;
        if parser.canonical_string()?.trim().is_empty() {
            return None;
        }
        parser.take(",\"residual_risk\":{\"likelihood\":")?;
        if !(1..=5).contains(&parser.canonical_u64()?) {
            return None;
        }
        parser.take(",\"impact\":")?;
        if !(1..=5).contains(&parser.canonical_u64()?) {
            return None;
        }
        parser.take("},\"review_gate\":")?;
        if !matches!(
            parser.canonical_string()?.as_str(),
            "E0a" | "E0b" | "E0c" | "E0d" | "E2" | "E4" | "E5" | "E6" | "E7"
        ) {
            return None;
        }
        parser.take("}")?;
    }
    parser.take("]}")?;
    (parser.is_finished() && triggers.len() == PROGRAM_RISK_REPORT_LOGICAL_ROWS)
        .then_some(HistoricRegisterV1 { triggers })
}

fn valid_historic_observed_unit(unit: &str, exact_bytes: u64) -> bool {
    let preview_bytes = u64::try_from(unit.len()).unwrap_or(u64::MAX);
    let preview_cap = u64::try_from(fs_govern::program_risks::MAX_OBSERVED_UNIT_PREVIEW_BYTES)
        .expect("program-risk preview cap fits u64");
    if exact_bytes <= preview_cap {
        exact_bytes == preview_bytes
    } else {
        let shortest_maximal_utf8_prefix = preview_cap.saturating_sub(3);
        (shortest_maximal_utf8_prefix..=preview_cap).contains(&preview_bytes)
            && exact_bytes > preview_bytes
    }
}

fn exact_historic_observed_unit(unit: &str, exact_bytes: u64) -> bool {
    u64::try_from(unit.len()).is_ok_and(|bytes| bytes == exact_bytes)
}

fn validate_historic_assessment_row_v1(
    parser: &mut ProgramRiskJsonCursor<'_>,
    expected_id: &str,
    expected_status: AssessmentStatus,
    trigger: &HistoricTriggerV1,
) -> Option<()> {
    parser.take("{\"id\":")?;
    if parser.canonical_string()? != expected_id {
        return None;
    }
    parser.take(",\"status\":")?;
    if parser.canonical_string()? != status_code_v1(expected_status) {
        return None;
    }
    parser.take(",\"observed_value\":")?;
    let observed_value = parser.canonical_optional_f64()?;
    parser.take(",\"observed_unit\":")?;
    let observed_unit = parser.canonical_optional_string()?;
    parser.take(",\"observed_unit_bytes\":")?;
    let observed_unit_bytes = parser.canonical_optional_u64()?;
    parser.take(",\"samples\":")?;
    let samples = parser.canonical_optional_u64()?;
    parser.take(",\"observation_count\":")?;
    let observation_count = parser.canonical_u64()?;
    parser.take(",\"trigger\":")?;
    if parse_historic_trigger_v1(parser)? != *trigger {
        return None;
    }
    parser.take("}")?;

    match expected_status {
        AssessmentStatus::Missing => {
            if observed_value.is_some()
                || observed_unit.is_some()
                || observed_unit_bytes.is_some()
                || samples.is_some()
                || observation_count != 0
            {
                return None;
            }
        }
        AssessmentStatus::Duplicate => {
            if observed_value.is_some()
                || observed_unit.is_some()
                || observed_unit_bytes.is_some()
                || samples.is_some()
                || observation_count < 2
            {
                return None;
            }
        }
        AssessmentStatus::NonFinite => {
            let (Some(unit), Some(unit_bytes), Some(_samples)) =
                (observed_unit.as_deref(), observed_unit_bytes, samples)
            else {
                return None;
            };
            if observed_value.is_some()
                || observation_count != 1
                || unit != trigger.unit
                || !exact_historic_observed_unit(unit, unit_bytes)
            {
                return None;
            }
        }
        AssessmentStatus::UnitMismatch => {
            let (Some(unit), Some(unit_bytes), Some(_samples)) =
                (observed_unit.as_deref(), observed_unit_bytes, samples)
            else {
                return None;
            };
            if observation_count != 1
                || (unit == trigger.unit && exact_historic_observed_unit(unit, unit_bytes))
                || !valid_historic_observed_unit(unit, unit_bytes)
            {
                return None;
            }
        }
        expected => {
            let (Some(value), Some(unit), Some(unit_bytes), Some(samples)) = (
                observed_value,
                observed_unit.as_deref(),
                observed_unit_bytes,
                samples,
            ) else {
                return None;
            };
            if observation_count != 1
                || !exact_historic_observed_unit(unit, unit_bytes)
                || trigger.assess(value, unit, samples) != expected
            {
                return None;
            }
        }
    }
    Some(())
}

fn validate_historic_report(
    bytes: &[u8],
    authority: fs_blake3::ContentHash,
    register: &HistoricRegisterV1,
    expected: &HistoricReportExpected<'_>,
) -> Result<(), SessionError> {
    let input = core::str::from_utf8(bytes)
        .map_err(|_| historic_corrupt(authority, "historic program-risk report is not UTF-8"))?;
    let mut parser = ProgramRiskJsonCursor { input, offset: 0 };
    parser
        .take("{\"schema\":")
        .ok_or_else(|| historic_corrupt(authority, "historic report lacks a canonical schema"))?;
    let schema = parser.canonical_string().ok_or_else(|| {
        historic_corrupt(
            authority,
            "historic report schema is not a canonical string",
        )
    })?;
    if schema != "frankensim.program-risk-session-report.v1" {
        return Err(historic_corrupt(
            authority,
            format!(
                "unsupported historic program-risk report schema {:?} ({} bytes)",
                bounded_utf8_prefix(&schema, 64),
                schema.len()
            ),
        ));
    }
    validate_historic_report_v1(&mut parser, register, expected).ok_or_else(|| {
        historic_corrupt(
            authority,
            "historic program-risk report v1 violates its frozen schema or terminal binding",
        )
    })
}

fn validate_historic_report_v1(
    parser: &mut ProgramRiskJsonCursor<'_>,
    register: &HistoricRegisterV1,
    expected: &HistoricReportExpected<'_>,
) -> Option<()> {
    parser.take(",\"session\":")?;
    if parser.canonical_u64()? != expected.session.0 {
        return None;
    }
    parser.take(",\"ledger_scope\":")?;
    if parser.canonical_string()? != expected.ledger_scope {
        return None;
    }
    parser.take(",\"session_open\":")?;
    if parser.canonical_string()? != expected.session_open.to_string() {
        return None;
    }
    parser.take(",\"register_artifact\":")?;
    if parser.canonical_string()? != expected.register_artifact.to_string() {
        return None;
    }
    parser.take(",\"logical_time\":")?;
    if parser.canonical_i64()? != expected.logical_time {
        return None;
    }
    parser.take(",\"generation\":")?;
    if parser.canonical_u64()? != expected.generation {
        return None;
    }
    parser.take(",\"assessment\":{\"schema\":")?;
    if parser.canonical_string()? != "frankensim.program-risk-assessment.v1" {
        return None;
    }
    let expected_alert_count = expected
        .statuses
        .iter()
        .filter(|status| status.is_alert())
        .count();
    parser.take(",\"all_clear\":")?;
    if parser.canonical_bool()? != (expected_alert_count == 0) {
        return None;
    }
    parser.take(",\"alert_count\":")?;
    if usize::try_from(parser.canonical_u64()?).ok()? != expected_alert_count {
        return None;
    }
    parser.take(",\"rows\":[")?;
    if register.triggers.len() != PROGRAM_RISK_REPORT_LOGICAL_ROWS {
        return None;
    }
    for (index, ((expected_id, expected_status), trigger)) in PROGRAM_RISK_V1_ROW_CODES
        .into_iter()
        .zip(expected.statuses.iter().copied())
        .zip(register.triggers.iter())
        .enumerate()
    {
        if index > 0 {
            parser.take(",")?;
        }
        validate_historic_assessment_row_v1(parser, expected_id, expected_status, trigger)?;
    }
    parser.take("]}}")?;
    parser.is_finished().then_some(())
}

#[allow(clippy::too_many_arguments)] // Exact historic byte lengths are part of the frozen lineage budget.
fn verify_lineage_shape(
    ledger: &Ledger,
    op: i64,
    report_id: ProgramRiskReportId,
    expected_ir: &str,
    register_artifact: fs_blake3::ContentHash,
    register_len: usize,
    report_artifact: fs_blake3::ContentHash,
    report_len: usize,
) -> Result<(), SessionError> {
    let stored = ledger
        .op(op)
        .map_err(|error| map_ledger("program-risk lineage op read", error))?
        .ok_or_else(|| SessionError::Persistence {
            what: format!("program-risk lineage op {op} is missing"),
        })?;
    let session_bytes = report_id.session.0.to_be_bytes();
    let expected_budget = lineage_budget(register_len, report_len);
    if stored.session.as_deref() != Some(session_bytes.as_slice())
        || stored.ir != expected_ir
        || stored.seed.as_slice() != report_id.content_hash.as_bytes()
        || stored.versions != PROGRAM_RISK_LINEAGE_VERSIONS_V1
        || stored.budget != expected_budget
        || stored.capability != PROGRAM_RISK_LINEAGE_CAPABILITY_V1
        || stored.t_start != 0
        || stored.t_end != Some(0)
        || stored.outcome.as_deref() != Some("ok")
        || stored.diag.is_some()
    {
        return Err(SessionError::Persistence {
            what: format!(
                "program-risk lineage op {op} does not match its exact session, IR, Five Explicits, deterministic timestamps, or completed outcome"
            ),
        });
    }

    let producers = ledger
        .artifact_producer_ops_bounded(&report_artifact, PROGRAM_RISK_REPORT_PRODUCER_CAP)
        .map_err(|error| map_ledger("program-risk report producer query", error))?;
    if producers.truncated || producers.op_ids.as_slice() != [op] {
        return Err(SessionError::Persistence {
            what: format!(
                "program-risk lineage op {op} is not the sole output producer of the report artifact"
            ),
        });
    }
    let context = ledger
        .op_execution_context(op)
        .map_err(|error| map_ledger("program-risk lineage execution context", error))?
        .ok_or_else(|| SessionError::Persistence {
            what: format!("program-risk lineage op {op} disappeared during verification"),
        })?;
    if context.exec_mode != ExecMode::Deterministic || context.branch != fs_ledger::MAIN_BRANCH {
        return Err(SessionError::Persistence {
            what: format!(
                "program-risk lineage op {op} was not recorded in deterministic mode on the main branch"
            ),
        });
    }

    let linked = ledger
        .op_artifact_edges_bounded(op, PROGRAM_RISK_LINEAGE_EDGE_CAP)
        .map_err(|error| map_ledger("program-risk exact lineage edge set", error))?;
    let expected_edges = [
        OpArtifactEdge {
            role: EdgeRole::In,
            artifact: register_artifact,
        },
        OpArtifactEdge {
            role: EdgeRole::Out,
            artifact: report_artifact,
        },
    ];
    if linked.truncated || linked.edges.as_slice() != expected_edges {
        return Err(SessionError::Persistence {
            what: format!(
                "program-risk lineage op {op} has extra, missing, duplicate, or role-mismatched artifact edges"
            ),
        });
    }
    Ok(())
}

fn verify_lineage_seals(
    ledger: &Ledger,
    op: i64,
    report_artifact: fs_blake3::ContentHash,
) -> Result<(), SessionError> {
    if lineage_seals_are_complete(ledger, op, report_artifact)? {
        return Ok(());
    }
    Err(SessionError::Persistence {
        what: format!(
            "program-risk lineage op {op} is missing one or both immutable lineage seals"
        ),
    })
}

fn lineage_seals_are_complete(
    ledger: &Ledger,
    op: i64,
    report_artifact: fs_blake3::ContentHash,
) -> Result<bool, SessionError> {
    let sealed = ledger
        .artifact_output_seal(&report_artifact)
        .map_err(|error| map_ledger("program-risk report output seal", error))?;
    if let Some(sealed_op) = sealed
        && sealed_op != op
    {
        return Err(SessionError::Persistence {
            what: format!(
                "program-risk report artifact is immutably sealed to lineage op {sealed_op}, not {op}"
            ),
        });
    }
    let edge_count = ledger
        .op_artifact_edge_seal(op)
        .map_err(|error| map_ledger("program-risk lineage edge-set seal", error))?;
    if let Some(sealed_count) = edge_count
        && sealed_count != PROGRAM_RISK_LINEAGE_EDGE_CAP
    {
        return Err(SessionError::Persistence {
            what: format!(
                "program-risk lineage op {op} is immutably sealed at {sealed_count} artifact edges, not the required exact two"
            ),
        });
    }
    Ok(sealed.is_some() && edge_count.is_some())
}

#[allow(clippy::too_many_arguments)] // Exact historic byte lengths are part of the frozen lineage budget.
fn verify_lineage_op(
    ledger: &Ledger,
    op: i64,
    report_id: ProgramRiskReportId,
    expected_ir: &str,
    register_artifact: fs_blake3::ContentHash,
    register_len: usize,
    report_artifact: fs_blake3::ContentHash,
    report_len: usize,
) -> Result<(), SessionError> {
    if lineage_seals_are_complete(ledger, op, report_artifact)? {
        verify_lineage_shape(
            ledger,
            op,
            report_id,
            expected_ir,
            register_artifact,
            register_len,
            report_artifact,
            report_len,
        )?;
        return Ok(());
    }
    if ledger.in_transaction() {
        return Err(SessionError::Persistence {
            what: "program-risk legacy lineage-seal adoption refuses a caller-owned transaction"
                .to_string(),
        });
    }

    ledger
        .begin()
        .map_err(|error| map_ledger("program-risk legacy lineage-seal adoption begin", error))?;
    let adoption = (|| -> Result<(), SessionError> {
        lineage_seals_are_complete(ledger, op, report_artifact)?;
        verify_lineage_shape(
            ledger,
            op,
            report_id,
            expected_ir,
            register_artifact,
            register_len,
            report_artifact,
            report_len,
        )?;
        ledger
            .seal_artifact_output(&report_artifact, op)
            .map_err(|error| map_ledger("program-risk legacy report output seal", error))?;
        ledger
            .seal_op_artifact_edges(op, PROGRAM_RISK_LINEAGE_EDGE_CAP)
            .map_err(|error| map_ledger("program-risk legacy lineage edge-set seal", error))?;
        Ok(())
    })();
    if let Err(error) = adoption {
        let _ = ledger.rollback();
        return Err(error);
    }
    if let Err(error) = ledger.commit() {
        let _ = ledger.rollback();
        return Err(map_ledger(
            "program-risk legacy lineage-seal adoption commit",
            error,
        ));
    }

    verify_lineage_shape(
        ledger,
        op,
        report_id,
        expected_ir,
        register_artifact,
        register_len,
        report_artifact,
        report_len,
    )?;
    verify_lineage_seals(ledger, op, report_artifact)
}

fn ensure_lineage_op(
    ledger: &Ledger,
    report_id: ProgramRiskReportId,
    register_artifact: fs_blake3::ContentHash,
    register_len: usize,
    report_artifact: fs_blake3::ContentHash,
    report_len: usize,
) -> Result<i64, SessionError> {
    let ir = lineage_ir(report_id, register_artifact, report_artifact);
    if ledger.in_transaction() {
        return Err(SessionError::Persistence {
            what: "program-risk lineage materialization refuses a caller-owned transaction"
                .to_string(),
        });
    }
    let budget = lineage_budget(register_len, report_len);
    let session_bytes = report_id.session.0.to_be_bytes();
    let explicits = FiveExplicits {
        seed: report_id.content_hash.as_bytes(),
        versions: PROGRAM_RISK_LINEAGE_VERSIONS_V1,
        budget: &budget,
        capability: PROGRAM_RISK_LINEAGE_CAPABILITY_V1,
    };
    ledger
        .begin()
        .map_err(|error| map_ledger("program-risk lineage begin", error))?;
    let write = (|| -> Result<i64, SessionError> {
        let producers = ledger
            .artifact_producer_ops_bounded(&report_artifact, PROGRAM_RISK_REPORT_PRODUCER_CAP)
            .map_err(|error| map_ledger("program-risk lineage lookup", error))?;
        if producers.truncated {
            return Err(SessionError::Persistence {
                what: "program-risk report artifact has multiple output producers; refusing ambiguous lineage"
                    .to_string(),
            });
        }
        if let Some(&existing) = producers.op_ids.first() {
            verify_lineage_shape(
                ledger,
                existing,
                report_id,
                &ir,
                register_artifact,
                register_len,
                report_artifact,
                report_len,
            )?;
            ledger
                .seal_artifact_output(&report_artifact, existing)
                .map_err(|error| map_ledger("program-risk report output seal", error))?;
            ledger
                .seal_op_artifact_edges(existing, PROGRAM_RISK_LINEAGE_EDGE_CAP)
                .map_err(|error| map_ledger("program-risk lineage edge-set seal", error))?;
            return Ok(existing);
        }

        let op = ledger
            .begin_op(Some(&session_bytes), &ir, &explicits, 0)
            .map_err(|error| map_ledger("program-risk lineage op begin", error))?;
        ledger
            .link(op, &register_artifact, EdgeRole::In)
            .map_err(|error| map_ledger("program-risk register input edge", error))?;
        ledger
            .link(op, &report_artifact, EdgeRole::Out)
            .map_err(|error| map_ledger("program-risk report output edge", error))?;
        ledger
            .seal_artifact_output(&report_artifact, op)
            .map_err(|error| map_ledger("program-risk report output seal", error))?;
        ledger
            .seal_op_artifact_edges(op, PROGRAM_RISK_LINEAGE_EDGE_CAP)
            .map_err(|error| map_ledger("program-risk lineage edge-set seal", error))?;
        ledger
            .finish_op(op, OpOutcome::Ok, None, 0)
            .map_err(|error| map_ledger("program-risk lineage finish", error))?;
        Ok(op)
    })();
    match write {
        Ok(op) => {
            if let Err(error) = ledger.commit() {
                let _ = ledger.rollback();
                return Err(map_ledger("program-risk lineage commit", error));
            }
            verify_lineage_op(
                ledger,
                op,
                report_id,
                &ir,
                register_artifact,
                register_len,
                report_artifact,
                report_len,
            )?;
            Ok(op)
        }
        Err(error) => {
            let _ = ledger.rollback();
            Err(error)
        }
    }
}

fn event_payload(receipt: &ProgramRiskReportReceipt) -> String {
    let alerts = receipt
        .alerts()
        .into_iter()
        .map(|alert| {
            format!(
                "{{\"id\":\"{}\",\"status\":\"{}\"}}",
                alert.id.code(),
                status_code_v1(alert.status)
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"schema\":\"frankensim.program-risk-report-event.v1\",\"report_id\":\"{}\",\"session_open\":\"{}\",\"register_artifact\":\"{}\",\"report_artifact\":\"{}\",\"lineage_op\":{},\"generation\":{},\"alert_count\":{},\"alerts\":[{}]}}",
        receipt.report_id.content_hash,
        receipt.report_id.session_open,
        receipt.register_artifact,
        receipt.report_artifact,
        receipt.lineage_op,
        receipt.generation,
        receipt.alert_count(),
        alerts,
    )
}

fn decode_existing_report(
    terminal: &fs_ledger::session_registry::StoredSessionTerminal,
    report_id: ProgramRiskReportId,
) -> Result<ProgramRiskReportReceipt, SessionError> {
    let payload = decode_terminal_receipt(report_id, &terminal.claim.payload)?;
    let receipt = decode_terminal_receipt(report_id, &terminal.receipt)?;
    if payload != receipt {
        return Err(SessionError::TerminalCorrupt {
            kind: crate::governor::recovery::KIND_PROGRAM_RISK_REPORT,
            authority: report_id.content_hash,
            detail: "program-risk claim payload and terminal receipt disagree".to_string(),
        });
    }
    Ok(receipt)
}

fn validate_existing_report(
    receipt: &ProgramRiskReportReceipt,
    register_artifact: fs_blake3::ContentHash,
    report_artifact: fs_blake3::ContentHash,
    logical_time: i64,
    statuses: [AssessmentStatus; PROGRAM_RISK_REPORT_LOGICAL_ROWS],
) -> Result<(), SessionError> {
    if receipt.register_artifact != register_artifact
        || receipt.report_artifact != report_artifact
        || receipt.logical_time != logical_time
        || receipt.statuses != statuses
    {
        return Err(SessionError::MutationConflict {
            kind: crate::governor::recovery::KIND_PROGRAM_RISK_REPORT,
            id: receipt.report_id.content_hash,
        });
    }
    Ok(())
}

impl Governor {
    /// Persist one drained session's PR-001--PR-012 snapshot.
    ///
    /// The session's ordinary scope must already have been flushed to this
    /// ledger, which proves the immutable open receipt and binds the scope to
    /// one physical sink. This method reports a quiescent point-in-time
    /// snapshot; the caller owns session-end coordination. It does not close,
    /// revoke, or flush the session token's ordinary mutation queues.
    ///
    /// Exact retry with the same logical time and observations appends zero
    /// terminal/event rows, including after the live session advances to a
    /// later execution generation once the session is quiescent again. The
    /// first schema-v9 retry of valid schema-v8 lineage may atomically install
    /// its two missing immutable seal rows. Any changed retry conflicts under
    /// the same singleton authority.
    ///
    /// # Errors
    /// Foreign/open-mismatched receipts, pending submissions or pauses,
    /// unflushed/wrong-sink scopes, malformed evidence, caller-owned ledger
    /// transactions, lineage failures, and persistence conflicts fail closed.
    #[allow(clippy::too_many_lines)] // Two-phase artifact/lineage then typed terminal publication.
    pub fn write_program_risk_session_end_report(
        &self,
        ledger: &Ledger,
        open_receipt: &SessionOpenReceipt,
        logical_time: i64,
        observations: &[ProgramRiskObservation<'_>],
    ) -> Result<ProgramRiskReportWrite, SessionError> {
        if ledger.in_transaction() {
            return Err(SessionError::Persistence {
                what: "program-risk report refuses a caller-owned ledger transaction".to_string(),
            });
        }
        let sink = ledger
            .checked_instance_id()
            .map_err(|error| map_ledger("program-risk ledger identity", error))?;
        self.ensure_program_risk_mutation_ready()?;
        let (session, ledger_scope, current_generation) =
            self.program_risk_session_context(open_receipt, sink)?;
        let report_id = report_id(self, open_receipt);
        if !program_risk_report_identity_transport_is_current(report_id) {
            return Err(SessionError::Persistence {
                what: "program-risk report authority failed its current identity transport"
                    .to_string(),
            });
        }
        let _publication = self.reserve_program_risk_publication(report_id.content_hash)?;

        let assessment = assess_program_risks(observations);
        let statuses = assessment_statuses(&assessment)?;
        let register_json = program_risk_register_json();
        let register_artifact = fs_blake3::hash_bytes(register_json.as_bytes());
        let existing_terminal = ledger
            .session_terminal(&report_id.content_hash)
            .map_err(|error| map_ledger("program-risk terminal lookup", error))?;
        let existing_receipt = existing_terminal
            .as_ref()
            .map(|terminal| decode_existing_report(terminal, report_id))
            .transpose()?;
        let report_generation = existing_receipt
            .as_ref()
            .map_or(current_generation, |receipt| receipt.generation);
        let report_json = report_artifact_json(
            session,
            &ledger_scope,
            open_receipt,
            register_artifact,
            logical_time,
            report_generation,
            &assessment,
        );
        let report_artifact = fs_blake3::hash_bytes(report_json.as_bytes());

        let receipt = if let Some(receipt) = existing_receipt {
            validate_existing_report(
                &receipt,
                register_artifact,
                report_artifact,
                logical_time,
                statuses,
            )?;
            verify_artifact(
                ledger,
                register_artifact,
                PROGRAM_RISK_REGISTER_ARTIFACT_KIND,
                register_json.as_bytes(),
            )?;
            verify_artifact(
                ledger,
                report_artifact,
                PROGRAM_RISK_SESSION_REPORT_ARTIFACT_KIND,
                report_json.as_bytes(),
            )?;
            let expected_ir = lineage_ir(report_id, register_artifact, report_artifact);
            verify_lineage_op(
                ledger,
                receipt.lineage_op,
                report_id,
                &expected_ir,
                register_artifact,
                register_json.len(),
                report_artifact,
                report_json.len(),
            )?;
            receipt
        } else {
            if ledger
                .session_mutation_claim(&report_id.content_hash)
                .map_err(|error| map_ledger("program-risk pending-claim lookup", error))?
                .is_some()
            {
                return Err(SessionError::IndeterminateMutation {
                    kind: crate::governor::recovery::KIND_PROGRAM_RISK_REPORT,
                    authority: report_id.content_hash,
                });
            }
            let register_put = ledger
                .put_artifact(
                    PROGRAM_RISK_REGISTER_ARTIFACT_KIND,
                    register_json.as_bytes(),
                    Some(PROGRAM_RISK_REGISTER_META),
                )
                .map_err(|error| map_ledger("program-risk register artifact write", error))?;
            let report_put = ledger
                .put_artifact(
                    PROGRAM_RISK_SESSION_REPORT_ARTIFACT_KIND,
                    report_json.as_bytes(),
                    Some(PROGRAM_RISK_SESSION_REPORT_META),
                )
                .map_err(|error| map_ledger("program-risk session artifact write", error))?;
            if register_put.hash != register_artifact || report_put.hash != report_artifact {
                return Err(SessionError::Persistence {
                    what: "program-risk artifact writer returned a non-canonical content hash"
                        .to_string(),
                });
            }
            let lineage_op = ensure_lineage_op(
                ledger,
                report_id,
                register_artifact,
                register_json.len(),
                report_artifact,
                report_json.len(),
            )?;
            ProgramRiskReportReceipt {
                report_id,
                register_artifact,
                report_artifact,
                lineage_op,
                logical_time,
                generation: report_generation,
                statuses,
            }
        };

        let terminal_bytes = encode_terminal_receipt(&receipt);
        let payload = event_payload(&receipt);
        let session_bytes = session.0.to_be_bytes();
        let event = EventRow {
            session: Some(&session_bytes),
            t: logical_time,
            kind: PROGRAM_RISK_REPORT_EVENT_KIND,
            payload: Some(&payload),
        };
        let events = [event];
        let claim = SessionMutationClaim {
            authority: report_id.content_hash,
            ledger_instance_id: sink,
            governor_hash: self.identity(),
            session_open_hash: open_receipt.content_hash(),
            kind: crate::governor::recovery::KIND_PROGRAM_RISK_REPORT,
            session: session.0,
            ledger_scope: &ledger_scope,
            generation: receipt.generation,
            causal_ordinal: None,
            payload: &terminal_bytes,
        };
        let group = SessionTerminalGroup {
            terminal: SessionTerminalRow {
                claim,
                permit: None,
                receipt: &terminal_bytes,
            },
            events: &events,
        };
        let groups = [group];
        let result = ledger
            .append_session_terminal_batch(&SessionTerminalBatch { groups: &groups })
            .map_err(|error| map_ledger("program-risk terminal batch", error))?;
        let disposition = match result {
            SessionTerminalBatchResult::Committed {
                terminals_inserted: 1,
                events_appended: 1,
                ..
            } => ProgramRiskReportDisposition::Committed,
            SessionTerminalBatchResult::Replayed { .. } => ProgramRiskReportDisposition::Replayed,
            other => {
                return Err(SessionError::Persistence {
                    what: format!(
                        "program-risk terminal batch returned unexpected cardinality {other:?}"
                    ),
                });
            }
        };
        Ok(ProgramRiskReportWrite {
            receipt,
            disposition,
        })
    }

    /// Recover one already-persisted program-risk report into a durable
    /// governor's global recovery fence.
    ///
    /// Recover the session open (and any lifecycle generations preceding this
    /// report) first. The method verifies the strict payload/receipt codec,
    /// complete terminal claim, owned event, both artifact bytes, and both
    /// lineage edges before marking the durable claim recovered.
    pub fn recover_program_risk_report(
        &self,
        ledger: &Ledger,
        open_receipt: &SessionOpenReceipt,
    ) -> Result<ProgramRiskReportReceipt, SessionError> {
        let sink = self.recovery_ledger(ledger)?;
        let (session, ledger_scope, current_generation) =
            self.program_risk_session_context(open_receipt, sink)?;
        let report_id = report_id(self, open_receipt);
        let _publication = self.reserve_program_risk_publication(report_id.content_hash)?;
        let terminal = ledger
            .session_terminal(&report_id.content_hash)
            .map_err(|error| map_ledger("recover program-risk terminal", error))?
            .ok_or(SessionError::RecoveryRequired {
                kind: crate::governor::recovery::KIND_PROGRAM_RISK_REPORT,
                authority: report_id.content_hash,
            })?;
        let receipt = decode_existing_report(&terminal, report_id)?;
        if receipt.generation > current_generation {
            return Err(SessionError::ProgramRiskReportGenerationAhead {
                id: session.0,
                report_generation: receipt.generation,
                recovered_generation: current_generation,
            });
        }
        let register_json = read_historic_artifact(
            ledger,
            receipt.register_artifact,
            PROGRAM_RISK_REGISTER_ARTIFACT_KIND,
        )?;
        let report_json = read_historic_artifact(
            ledger,
            receipt.report_artifact,
            PROGRAM_RISK_SESSION_REPORT_ARTIFACT_KIND,
        )?;
        let historic_register =
            decode_historic_register(&register_json, receipt.register_artifact)?;
        validate_historic_report(
            &report_json,
            receipt.report_artifact,
            &historic_register,
            &HistoricReportExpected {
                session,
                ledger_scope: &ledger_scope,
                session_open: open_receipt.content_hash(),
                register_artifact: receipt.register_artifact,
                logical_time: receipt.logical_time,
                generation: receipt.generation,
                statuses: &receipt.statuses,
            },
        )?;
        let expected_ir = lineage_ir(
            report_id,
            receipt.register_artifact,
            receipt.report_artifact,
        );
        verify_lineage_op(
            ledger,
            receipt.lineage_op,
            report_id,
            &expected_ir,
            receipt.register_artifact,
            register_json.len(),
            receipt.report_artifact,
            report_json.len(),
        )?;
        let event_payload = event_payload(&receipt);
        let session_bytes = session.0.to_be_bytes();
        let expected_events = [EventRow {
            session: Some(&session_bytes),
            t: receipt.logical_time,
            kind: PROGRAM_RISK_REPORT_EVENT_KIND,
            payload: Some(&event_payload),
        }];
        let terminal_bytes = encode_terminal_receipt(&receipt);
        self.validate_recovered_terminal(
            &terminal,
            sink,
            crate::governor::recovery::KIND_PROGRAM_RISK_REPORT,
            report_id.content_hash,
            session,
            &ledger_scope,
            open_receipt.content_hash(),
            receipt.generation,
            None,
            &terminal_bytes,
            &expected_events,
        )?;
        self.mark_program_risk_report_recovered(report_id.content_hash);
        Ok(receipt)
    }
}

#[cfg(test)]
mod identity_tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};

    use crate::{CapabilityToken, DurableGovernorNonce};

    fn hash(byte: u8) -> fs_blake3::ContentHash {
        fs_blake3::ContentHash([byte; 32])
    }

    fn legacy_token(session: SessionId, ledger_scope: &str) -> CapabilityToken {
        CapabilityToken {
            session,
            ops: vec!["governance.*".to_string()],
            core_s: 60.0,
            mem_bytes: 1024 * 1024,
            wall_s: 60.0,
            cores: 1,
            ledger_scope: ledger_scope.to_string(),
        }
    }

    fn durable_legacy_path(case: &str) -> String {
        static NEXT: AtomicU64 = AtomicU64::new(1);
        let ordinal = NEXT.fetch_add(1, Ordering::Relaxed);
        std::env::temp_dir()
            .join(format!(
                "fs-session-program-risk-v8-migration-{}-{ordinal}-{case}.ledger",
                std::process::id()
            ))
            .to_string_lossy()
            .into_owned()
    }

    fn non_seal_row_counts(ledger: &Ledger) -> [u64; 6] {
        [
            ledger.table_count("artifacts").expect("artifact count"),
            ledger.table_count("ops").expect("op count"),
            ledger.table_count("edges").expect("edge count"),
            ledger
                .table_count("session_terminals")
                .expect("terminal count"),
            ledger
                .table_count("session_terminal_events")
                .expect("terminal-event count"),
            ledger.table_count("events").expect("event count"),
        ]
    }

    fn durable_report_row_counts(ledger: &Ledger) -> [u64; 8] {
        let [artifacts, ops, edges, terminals, terminal_events, events] =
            non_seal_row_counts(ledger);
        [
            artifacts,
            ops,
            edges,
            terminals,
            terminal_events,
            events,
            ledger
                .table_count("artifact_output_seals")
                .expect("output-seal count"),
            ledger
                .table_count("op_artifact_edge_seals")
                .expect("edge-set-seal count"),
        ]
    }

    /// Reproduce the populated row state immediately after a schema-v8 ledger
    /// migrates to v9: the PPVS artifacts, exact lineage, terminal, and event
    /// exist, while the newly-created seal tables are still empty.
    #[allow(clippy::too_many_lines)] // One exact populated legacy terminal fixture.
    fn persist_migrated_v8_style_report(
        governor: &Governor,
        ledger: &Ledger,
        open_receipt: &SessionOpenReceipt,
        logical_time: i64,
        observations: &[ProgramRiskObservation<'_>],
    ) -> ProgramRiskReportReceipt {
        let sink = ledger
            .checked_instance_id()
            .expect("legacy ledger identity");
        let (session, ledger_scope, generation) = governor
            .program_risk_session_context(open_receipt, sink)
            .expect("legacy report context");
        let report_id = report_id(governor, open_receipt);
        let assessment = assess_program_risks(observations);
        let statuses = assessment_statuses(&assessment).expect("twelve legacy statuses");
        let register_json = program_risk_register_json();
        let register_artifact = fs_blake3::hash_bytes(register_json.as_bytes());
        let report_json = report_artifact_json(
            session,
            &ledger_scope,
            open_receipt,
            register_artifact,
            logical_time,
            generation,
            &assessment,
        );
        let report_artifact = fs_blake3::hash_bytes(report_json.as_bytes());
        assert_eq!(
            ledger
                .put_artifact(
                    PROGRAM_RISK_REGISTER_ARTIFACT_KIND,
                    register_json.as_bytes(),
                    Some(PROGRAM_RISK_REGISTER_META),
                )
                .expect("legacy register artifact")
                .hash,
            register_artifact
        );
        assert_eq!(
            ledger
                .put_artifact(
                    PROGRAM_RISK_SESSION_REPORT_ARTIFACT_KIND,
                    report_json.as_bytes(),
                    Some(PROGRAM_RISK_SESSION_REPORT_META),
                )
                .expect("legacy report artifact")
                .hash,
            report_artifact
        );

        let ir = lineage_ir(report_id, register_artifact, report_artifact);
        let budget = lineage_budget(register_json.len(), report_json.len());
        let session_bytes = session.0.to_be_bytes();
        let explicits = FiveExplicits {
            seed: report_id.content_hash.as_bytes(),
            versions: PROGRAM_RISK_LINEAGE_VERSIONS_V1,
            budget: &budget,
            capability: PROGRAM_RISK_LINEAGE_CAPABILITY_V1,
        };
        let lineage_op = ledger
            .begin_op(Some(&session_bytes), &ir, &explicits, 0)
            .expect("legacy lineage op");
        ledger
            .link(lineage_op, &register_artifact, EdgeRole::In)
            .expect("legacy register edge");
        ledger
            .link(lineage_op, &report_artifact, EdgeRole::Out)
            .expect("legacy report edge");
        ledger
            .finish_op(lineage_op, OpOutcome::Ok, None, 0)
            .expect("legacy lineage finish");

        let receipt = ProgramRiskReportReceipt {
            report_id,
            register_artifact,
            report_artifact,
            lineage_op,
            logical_time,
            generation,
            statuses,
        };
        let terminal_bytes = encode_terminal_receipt(&receipt);
        let payload = event_payload(&receipt);
        let event = EventRow {
            session: Some(&session_bytes),
            t: logical_time,
            kind: PROGRAM_RISK_REPORT_EVENT_KIND,
            payload: Some(&payload),
        };
        let events = [event];
        let claim = SessionMutationClaim {
            authority: report_id.content_hash,
            ledger_instance_id: sink,
            governor_hash: governor.identity(),
            session_open_hash: open_receipt.content_hash(),
            kind: crate::governor::recovery::KIND_PROGRAM_RISK_REPORT,
            session: session.0,
            ledger_scope: &ledger_scope,
            generation,
            causal_ordinal: None,
            payload: &terminal_bytes,
        };
        let group = SessionTerminalGroup {
            terminal: SessionTerminalRow {
                claim,
                permit: None,
                receipt: &terminal_bytes,
            },
            events: &events,
        };
        let groups = [group];
        assert!(matches!(
            ledger
                .append_session_terminal_batch(&SessionTerminalBatch { groups: &groups })
                .expect("legacy terminal batch"),
            SessionTerminalBatchResult::Committed {
                terminals_inserted: 1,
                events_appended: 1,
                ..
            }
        ));
        assert_eq!(
            ledger
                .artifact_output_seal(&report_artifact)
                .expect("legacy output seal absence"),
            None
        );
        assert_eq!(
            ledger
                .op_artifact_edge_seal(lineage_op)
                .expect("legacy edge-set seal absence"),
            None
        );
        receipt
    }

    #[test]
    fn migrated_v8_populated_report_replay_adopts_only_missing_v9_seals() {
        let ledger = Ledger::open(":memory:").expect("migrated-v8 replay ledger");
        let governor = Governor::new();
        let session = SessionId(91_001);
        let open_id = governor
            .session_open_id(session, "migrated-v8-replay-open")
            .expect("open authority");
        let open = governor
            .open_session_declared(open_id, legacy_token(session, "migrated-v8-replay"))
            .expect("open session");
        governor
            .flush_scope_to_ledger(&open.flush_permit(), &ledger)
            .expect("persist open prerequisite");
        let observations: [ProgramRiskObservation<'static>; 0] = [];
        let legacy = persist_migrated_v8_style_report(&governor, &ledger, &open, 51, &observations);
        let counts = non_seal_row_counts(&ledger);

        let replay = governor
            .write_program_risk_session_end_report(&ledger, &open, 51, &observations)
            .expect("first schema-v9 replay adopts seals");
        assert_eq!(replay.disposition, ProgramRiskReportDisposition::Replayed);
        assert_eq!(replay.receipt, legacy);
        assert_eq!(non_seal_row_counts(&ledger), counts);
        assert_eq!(
            ledger
                .artifact_output_seal(&legacy.report_artifact)
                .expect("adopted output seal"),
            Some(legacy.lineage_op)
        );
        assert_eq!(
            ledger
                .op_artifact_edge_seal(legacy.lineage_op)
                .expect("adopted edge-set seal"),
            Some(PROGRAM_RISK_LINEAGE_EDGE_CAP)
        );
    }

    #[test]
    fn migrated_v8_populated_report_recovery_adopts_seals_before_acceptance() {
        let path = durable_legacy_path("recovery");
        let nonce = DurableGovernorNonce::from_bytes([0xD8; 32]);
        let session = SessionId(91_002);
        let capability = legacy_token(session, "migrated-v8-recovery");
        let observations: [ProgramRiskObservation<'static>; 0] = [];

        let ledger = Ledger::open(&path).expect("legacy durable ledger");
        let governor = Governor::new_durable(&ledger, nonce).expect("legacy durable governor");
        let open_id = governor
            .session_open_id(session, "migrated-v8-recovery-open")
            .expect("open authority");
        let open = governor
            .open_session_declared(open_id, capability.clone())
            .expect("open durable session");
        governor
            .flush_scope_to_ledger(&open.flush_permit(), &ledger)
            .expect("persist durable open");
        let legacy = persist_migrated_v8_style_report(&governor, &ledger, &open, 57, &observations);
        let counts = non_seal_row_counts(&ledger);
        drop(governor);
        drop(ledger);

        let ledger = Ledger::open(&path).expect("reopened migrated-v8 ledger");
        let governor = Governor::new_durable(&ledger, nonce).expect("reopened governor");
        let recovered_open = governor
            .recover_open(&ledger, open_id, capability, None)
            .expect("recover durable open first");
        let recovered = governor
            .recover_program_risk_report(&ledger, &recovered_open)
            .expect("schema-v9 recovery adopts legacy seals");
        assert_eq!(recovered, legacy);
        assert_eq!(non_seal_row_counts(&ledger), counts);
        assert_eq!(
            ledger
                .artifact_output_seal(&legacy.report_artifact)
                .expect("recovery output seal"),
            Some(legacy.lineage_op)
        );
        assert_eq!(
            ledger
                .op_artifact_edge_seal(legacy.lineage_op)
                .expect("recovery edge-set seal"),
            Some(PROGRAM_RISK_LINEAGE_EDGE_CAP)
        );

        let sealed_counts = durable_report_row_counts(&ledger);
        let replay = governor
            .write_program_risk_session_end_report(&ledger, &recovered_open, 57, &observations)
            .expect("post-recovery exact replay");
        assert_eq!(replay.disposition, ProgramRiskReportDisposition::Replayed);
        assert_eq!(replay.receipt, legacy);
        assert_eq!(durable_report_row_counts(&ledger), sealed_counts);
    }

    #[test]
    fn program_risk_report_identity_fields_move_independently() {
        let base = ProgramRiskReportIdIdentitySource {
            governor_id: hash(1),
            session: SessionId(2),
            session_open: hash(3),
        };
        let identity = program_risk_report_authority_identity(&base);
        for changed in [
            ProgramRiskReportIdIdentitySource {
                governor_id: hash(4),
                ..base
            },
            ProgramRiskReportIdIdentitySource {
                session: SessionId(5),
                ..base
            },
            ProgramRiskReportIdIdentitySource {
                session_open: hash(6),
                ..base
            },
        ] {
            assert_ne!(identity, program_risk_report_authority_identity(&changed));
        }

        let mut canonical = Vec::with_capacity(73);
        canonical.extend_from_slice(base.governor_id.as_bytes());
        canonical.extend_from_slice(&base.session.0.to_le_bytes());
        canonical.extend_from_slice(base.session_open.as_bytes());
        canonical.push(PROGRAM_RISK_REPORT_SLOT_TAG);
        assert_ne!(
            identity,
            fs_blake3::hash_domain(
                "org.frankensim.fs-session.program-risk-report-id.alternate.v2",
                &canonical,
            ),
            "digest domain is semantic",
        );
        assert_ne!(
            identity,
            fs_blake3::hash_domain(
                "org.frankensim.fs-session.program-risk-report-id.v3",
                &canonical,
            ),
            "identity version is semantic",
        );
        let mut reordered = Vec::with_capacity(73);
        reordered.extend_from_slice(base.session_open.as_bytes());
        reordered.extend_from_slice(&base.session.0.to_le_bytes());
        reordered.extend_from_slice(base.governor_id.as_bytes());
        reordered.push(PROGRAM_RISK_REPORT_SLOT_TAG);
        assert_ne!(
            identity,
            fs_blake3::hash_domain(PROGRAM_RISK_REPORT_ID_DOMAIN, &reordered),
            "canonical field order is semantic",
        );
        canonical[72] = PROGRAM_RISK_REPORT_SLOT_TAG + 1;
        assert_ne!(
            identity,
            fs_blake3::hash_domain(PROGRAM_RISK_REPORT_ID_DOMAIN, &canonical),
            "singleton slot is semantic",
        );
    }

    #[test]
    fn program_risk_report_identity_version_and_transport_fail_closed() {
        assert_eq!(PROGRAM_RISK_REPORT_IDENTITY_VERSION, 2);
        assert_eq!(PROGRAM_RISK_REPORT_SLOT_TAG, 1);
        let source = ProgramRiskReportIdIdentitySource {
            governor_id: hash(1),
            session: SessionId(2),
            session_open: hash(3),
        };
        let id = ProgramRiskReportId {
            governor_id: source.governor_id,
            session: source.session,
            session_open: source.session_open,
            content_hash: program_risk_report_authority_identity(&source),
        };
        assert!(program_risk_report_identity_transport_is_current(id));
        let altered = ProgramRiskReportId {
            content_hash: hash(9),
            ..id
        };
        assert!(!program_risk_report_identity_transport_is_current(altered));
    }

    #[test]
    fn program_risk_terminal_codec_round_trips_generation_and_statuses() {
        let source = ProgramRiskReportIdIdentitySource {
            governor_id: hash(1),
            session: SessionId(2),
            session_open: hash(3),
        };
        let report_id = ProgramRiskReportId {
            governor_id: source.governor_id,
            session: source.session,
            session_open: source.session_open,
            content_hash: program_risk_report_authority_identity(&source),
        };
        let statuses = [
            AssessmentStatus::Clear,
            AssessmentStatus::Triggered,
            AssessmentStatus::Missing,
            AssessmentStatus::Duplicate,
            AssessmentStatus::NonFinite,
            AssessmentStatus::UnitMismatch,
            AssessmentStatus::OutOfRange,
            AssessmentStatus::UnderSampled,
            AssessmentStatus::Clear,
            AssessmentStatus::Clear,
            AssessmentStatus::Clear,
            AssessmentStatus::Clear,
        ];
        let receipt = ProgramRiskReportReceipt {
            report_id,
            register_artifact: hash(4),
            report_artifact: hash(5),
            lineage_op: 6,
            logical_time: -7,
            generation: 8,
            statuses,
        };
        let encoded = encode_terminal_receipt(&receipt);
        assert_eq!(PROGRAM_RISK_REPORT_CODEC_VERSION, 1);
        assert_eq!(PROGRAM_RISK_REPORT_ROW_ORDER_VERSION, 1);
        assert_eq!(PROGRAM_RISK_REPORT_STATUS_TAG_VERSION, 1);
        assert_eq!(PROGRAM_RISK_REPORT_CODEC_V1_BYTES, 104);
        assert_eq!(encoded.len(), PROGRAM_RISK_REPORT_CODEC_V1_BYTES);
        assert_eq!(&encoded[0..4], &1_u32.to_le_bytes());
        assert_eq!(&encoded[4..36], hash(4).as_bytes());
        assert_eq!(&encoded[36..68], hash(5).as_bytes());
        assert_eq!(&encoded[68..76], &6_i64.to_le_bytes());
        assert_eq!(&encoded[76..84], &(-7_i64).to_le_bytes());
        assert_eq!(&encoded[84..92], &8_u64.to_le_bytes());
        assert_eq!(&encoded[92..100], &[0, 1, 2, 3, 4, 5, 6, 7]);
        assert_eq!(
            decode_terminal_receipt(report_id, &encoded).expect("current codec"),
            receipt
        );

        let mut future = encoded.clone();
        future[0..4].copy_from_slice(&2_u32.to_le_bytes());
        assert!(matches!(
            decode_terminal_receipt(report_id, &future),
            Err(SessionError::UnsupportedTerminalSchema {
                found: 2,
                supported: 1
            })
        ));

        let mut retired = encoded;
        retired[0..4].copy_from_slice(&0_u32.to_le_bytes());
        assert!(matches!(
            decode_terminal_receipt(report_id, &retired),
            Err(SessionError::TerminalCorrupt { .. })
        ));
    }

    #[test]
    fn program_risk_terminal_row_and_status_mappings_are_frozen() {
        assert_eq!(
            PROGRAM_RISK_V1_ROW_CODES,
            [
                "PR-001", "PR-002", "PR-003", "PR-004", "PR-005", "PR-006", "PR-007", "PR-008",
                "PR-009", "PR-010", "PR-011", "PR-012",
            ]
        );
        assert_eq!(
            ProgramRiskId::ALL.map(ProgramRiskId::code),
            PROGRAM_RISK_V1_ROW_CODES
        );
        assert_eq!(ProgramRiskId::ALL, PROGRAM_RISK_V1_IDS);
        for (status, tag) in [
            (AssessmentStatus::Clear, 0),
            (AssessmentStatus::Triggered, 1),
            (AssessmentStatus::Missing, 2),
            (AssessmentStatus::Duplicate, 3),
            (AssessmentStatus::NonFinite, 4),
            (AssessmentStatus::UnitMismatch, 5),
            (AssessmentStatus::OutOfRange, 6),
            (AssessmentStatus::UnderSampled, 7),
        ] {
            assert_eq!(status_tag(status), tag);
            assert_eq!(
                status_from_tag(tag, hash(9)).expect("v1 status tag"),
                status
            );
            assert_eq!(status_code_v1(status), status.code());
        }
    }

    #[test]
    fn program_risk_publication_reservation_is_process_local_and_scoped() {
        let governor = Governor::new();
        let authority = hash(9);
        let first = governor
            .reserve_program_risk_publication(authority)
            .expect("first reservation");
        assert!(matches!(
            governor.reserve_program_risk_publication(authority),
            Err(SessionError::ProgramRiskReportInFlight { id }) if id == authority
        ));
        drop(first);
        let _retry = governor
            .reserve_program_risk_publication(authority)
            .expect("reservation released on guard drop");
    }

    #[test]
    fn historic_v1_report_is_bound_to_terminal_statuses_and_register_triggers() {
        let register_json = program_risk_register_json();
        let register_artifact = fs_blake3::hash_bytes(register_json.as_bytes());
        let register = decode_historic_register(register_json.as_bytes(), register_artifact)
            .expect("current register satisfies frozen v1 decoder");
        let assessment = assess_program_risks(&[]);
        let statuses = assessment_statuses(&assessment).expect("twelve statuses");
        let session = SessionId(2);
        let ledger_scope = "historic\\scope";
        let session_open = hash(3);
        let report_json = format!(
            "{{\"schema\":\"frankensim.program-risk-session-report.v1\",\"session\":{},\"ledger_scope\":\"{}\",\"session_open\":\"{}\",\"register_artifact\":\"{}\",\"logical_time\":-7,\"generation\":8,\"assessment\":{}}}",
            session.0,
            json_escape(ledger_scope),
            session_open,
            register_artifact,
            assessment.to_json(),
        );
        let report_artifact = fs_blake3::hash_bytes(report_json.as_bytes());
        let expected = HistoricReportExpected {
            session,
            ledger_scope,
            session_open,
            register_artifact,
            logical_time: -7,
            generation: 8,
            statuses: &statuses,
        };
        validate_historic_report(
            report_json.as_bytes(),
            report_artifact,
            &register,
            &expected,
        )
        .expect("current report satisfies frozen v1 decoder");

        let altered = report_json.replacen("\"status\":\"missing\"", "\"status\":\"triggered\"", 1);
        assert!(matches!(
            validate_historic_report(
                altered.as_bytes(),
                fs_blake3::hash_bytes(altered.as_bytes()),
                &register,
                &expected,
            ),
            Err(SessionError::TerminalCorrupt { .. })
        ));

        let unit_mismatch_assessment = assess_program_risks(&[ProgramRiskObservation {
            id: ProgramRiskId::Pr001,
            value: f64::NAN,
            unit: "",
            samples: 5,
        }]);
        let unit_mismatch_statuses =
            assessment_statuses(&unit_mismatch_assessment).expect("twelve statuses");
        assert_eq!(
            unit_mismatch_statuses[0],
            AssessmentStatus::UnitMismatch,
            "unit mismatch takes precedence over non-finite evidence"
        );
        let unit_mismatch_report = format!(
            "{{\"schema\":\"frankensim.program-risk-session-report.v1\",\"session\":{},\"ledger_scope\":\"{}\",\"session_open\":\"{}\",\"register_artifact\":\"{}\",\"logical_time\":-7,\"generation\":8,\"assessment\":{}}}",
            session.0,
            json_escape(ledger_scope),
            session_open,
            register_artifact,
            unit_mismatch_assessment.to_json(),
        );
        let unit_mismatch_expected = HistoricReportExpected {
            statuses: &unit_mismatch_statuses,
            ..expected
        };
        validate_historic_report(
            unit_mismatch_report.as_bytes(),
            fs_blake3::hash_bytes(unit_mismatch_report.as_bytes()),
            &register,
            &unit_mismatch_expected,
        )
        .expect("non-finite unit mismatch remains recoverable under v1 precedence");
    }
}

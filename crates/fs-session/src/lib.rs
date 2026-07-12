//! fs-session (plan §11.3): sessions, capability tokens, and the resource
//! GOVERNOR — budgets are ENFORCED, not advisory — plus the agent-proofing
//! trio: idempotency keys (a retry cannot double-spend), `estimate()` dry
//! runs (plan before you spend), and errors as GUIDANCE ("a refusal that
//! teaches is worth ten silent successes").
//!
//! Layer: L6 (HELM). Threading contract: the governor's hot paths are
//! `Send + Sync` (in-memory, mutex-guarded) so enforcement and idempotency
//! survive concurrent submission storms; ledger persistence is an explicit
//! single-threaded `flush_scope_to_ledger` step because fsqlite connections are
//! `!Send` by design.

pub mod estimate;
pub mod gemm_tune;
pub mod governor;
pub mod guidance;
pub mod token;

pub use estimate::{
    CalibrationHealth, CalibrationPolicy, CalibrationReport, Estimate, ZeroPredictionSummary,
    estimate,
};
pub use gemm_tune::{
    GEMM_DEPGRAPH_RECEIPT_DOMAIN, GEMM_TUNE_METADATA_PLAN_SCHEMA, GEMM_TUNE_ROW_RECEIPT_DOMAIN,
    GEMM_TUNER_SCHEMA_VERSION, GemmDispatch, GemmExecutionReceipt, GemmGraphEvidenceClass,
    GemmMemoryReceipt, GemmPanelReceipt, GemmTuneBuildEvidence, GemmTuneCache, GemmTuneError,
    ValidatedGemmTuneRow, gemm_f64_session, gemm_f64_session_budgeted, gemm_f64_session_with_pool,
    gemm_f64_session_with_pool_budgeted, gemm_f64_session_with_pool_declared,
    gemm_f64_session_with_pool_declared_budgeted, gemm_kernel_key, gemm_shape_class,
    gemm_tune_build_evidence, gemm_tune_key, gemm_tune_key_budgeted, gemm_tune_key_with_pool,
    gemm_tune_key_with_pool_budgeted, gemm_tune_metadata_plan_bytes,
};
pub use governor::{
    Charge, DegradationEvent, DegradationStep, Enforcement, FlushReport, Governor,
    MAX_CHECKPOINT_CLAIM_BYTES, MAX_DEGRADATION_EVENTS_PER_SCOPE, MAX_EVENT_PAGE_ROWS,
    MAX_FLUSH_ENCODED_BYTES, MAX_FLUSH_ROWS, MAX_IDEMPOTENCY_INPUT_BYTES,
    MAX_IDEMPOTENCY_KEYS_PER_SESSION, MAX_RETAINED_BYTES_PER_GOVERNOR,
    MAX_RETAINED_BYTES_PER_SCOPE, MAX_RETAINED_EVIDENCE_BYTES, MAX_SESSIONS_PER_GOVERNOR,
    MAX_SESSIONS_PER_SCOPE, PauseAcknowledgement, PauseRequestId, RetainedEvidence,
    ScopeFlushPermit, StepPhase, SubmissionReceipt, SubmitOutcome,
};
pub use guidance::Guidance;
pub use token::{
    CapabilityToken, MAX_CAPABILITY_OP_BYTES, MAX_CAPABILITY_OPS, MAX_CAPABILITY_TOTAL_OP_BYTES,
    MAX_LEDGER_SCOPE_BYTES, SessionId,
};

use core::fmt;

/// Crate version (compile-time stamp).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Structured session failures (Decalogue P10).
#[derive(Debug, Clone, PartialEq)]
pub enum SessionError {
    /// The session id is unknown to the governor.
    UnknownSession {
        /// The id.
        id: u64,
    },
    /// A session id was registered more than once. Session identity is
    /// immutable: replacing a token would let new authority inherit old
    /// meters, pause state, and idempotency receipts.
    SessionAlreadyOpen {
        /// The duplicate id.
        id: u64,
    },
    /// A ledger scope was not a canonical bounded authority string.
    InvalidLedgerScope {
        /// UTF-8-safe prefix of the rejected string, bounded to the maximum
        /// admitted scope length.
        scope_preview: String,
        /// Exact byte length of the rejected string.
        scope_bytes: usize,
        /// Canonical scope grammar.
        requirement: &'static str,
    },
    /// One operator grant is not a bounded canonical authority string.
    InvalidOperatorGrant {
        /// Position in the token's operator list.
        index: usize,
        /// Bounded diagnostic prefix.
        grant_preview: String,
        /// Exact input byte length.
        grant_bytes: usize,
        /// Canonical grant grammar.
        requirement: &'static str,
    },
    /// A token repeats one operator grant.
    DuplicateOperatorGrant {
        /// Exact already-bounded duplicate.
        grant: String,
    },
    /// No open session carries the requested ledger scope.
    UnknownLedgerScope {
        /// Requested exact scope.
        scope: String,
    },
    /// A scope already persisted to a different ledger sink.
    LedgerScopeSinkMismatch {
        /// Scope whose history would be split.
        scope: String,
        /// Sink bound by the first successful non-empty flush.
        bound_sink: fs_ledger::LedgerInstanceId,
        /// Rejected sink.
        attempted_sink: fs_ledger::LedgerInstanceId,
    },
    /// A scoped flush permit was minted by a different governor.
    ScopePermitMismatch {
        /// Exact bounded scope carried by the foreign permit.
        scope: String,
    },
    /// A flush for this scope is already outside the state lock performing
    /// ledger I/O; another flush must retry rather than race its cursors.
    ScopeFlushInFlight {
        /// Exact bounded scope.
        scope: String,
    },
    /// A deterministic governor collection, payload, or ordinal bound was
    /// reached. Refusal happens before partial state mutation.
    LimitExceeded {
        /// Bounded resource name.
        resource: &'static str,
        /// Configured maximum.
        limit: usize,
        /// Exact observation or conservative lower bound.
        observed_at_least: usize,
    },
    /// A resource grant, charge, or accumulated meter is outside its valid
    /// finite, non-negative domain.
    InvalidResource {
        /// The resource field.
        resource: &'static str,
        /// The rejected value.
        value: f64,
        /// The required domain.
        requirement: &'static str,
    },
    /// A submission failed structurally (parse/admission).
    Submission {
        /// Diagnosis.
        what: String,
    },
    /// Ledger persistence failed.
    Persistence {
        /// Diagnosis.
        what: String,
    },
    /// A memory-pressure level outside the declared ladder 1..=3
    /// (bead gp3.13: out-of-ladder levels are refused, never clamped).
    InvalidPressureLevel {
        /// The rejected level.
        level: u8,
    },
    /// Level-3 pressure targeted a session opened without a bound
    /// cancellation gate — a pause that cannot reach the computation
    /// is refused, not ledgered (bead gp3.13).
    UngatedSession {
        /// The id.
        id: u64,
    },
    /// A caller attempted to bind an already-requested cancellation gate.
    PreRequestedGate {
        /// The session that would have inherited stale cancellation.
        id: u64,
    },
    /// An acknowledgement request belongs to another governor or does not
    /// match the session's pending/completed generation.
    PauseRequestMismatch {
        /// Session named by the opaque request.
        id: u64,
        /// Request ordinal carried by the stale/foreign authority.
        requested_ordinal: i64,
    },
    /// A completed request was replayed with different checkpoint evidence.
    PauseAcknowledgementConflict {
        /// Session whose terminal acknowledgement cannot be replaced.
        id: u64,
        /// Completed request ordinal.
        requested_ordinal: i64,
    },
    /// Pressure arrived before a fresh resume gate was explicitly activated.
    ResumeNotActivated {
        /// Session awaiting activation.
        id: u64,
        /// Fresh gate generation awaiting activation.
        generation: u64,
    },
    /// A supplied acknowledgement is stale, altered, or from another governor.
    ResumeAcknowledgementMismatch {
        /// Session named by the acknowledgement.
        id: u64,
    },
    /// The fresh resume gate was requested before activation completed.
    ResumeGateAlreadyRequested {
        /// Session whose gate is already cancelled.
        id: u64,
        /// Affected gate generation.
        generation: u64,
    },
    /// A pressure transition was requested while an earlier pause request
    /// still awaits its checkpoint acknowledgement.
    PauseAlreadyPending {
        /// The session id.
        id: u64,
        /// Ordinal of the still-pending request.
        requested_ordinal: i64,
    },
}

impl fmt::Display for SessionError {
    #[allow(clippy::too_many_lines)] // Exhaustive rendering stays adjacent to the error variants.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SessionError::UnknownSession { id } => write!(f, "unknown session {id}"),
            SessionError::SessionAlreadyOpen { id } => write!(
                f,
                "session {id} is already open; capability tokens are immutable and the existing \
                 session state was left unchanged"
            ),
            SessionError::InvalidLedgerScope {
                scope_preview,
                scope_bytes,
                requirement,
            } => write!(
                f,
                "invalid ledger scope {scope_preview:?} (input bytes: {scope_bytes}): {requirement}; session and flush state were not mutated"
            ),
            SessionError::InvalidOperatorGrant {
                index,
                grant_preview,
                grant_bytes,
                requirement,
            } => write!(
                f,
                "invalid operator grant {index} {grant_preview:?} (input bytes: {grant_bytes}): {requirement}; session authority was not registered"
            ),
            SessionError::DuplicateOperatorGrant { grant } => write!(
                f,
                "duplicate operator grant {grant:?}; session authority was not registered"
            ),
            SessionError::UnknownLedgerScope { scope } => write!(
                f,
                "unknown ledger scope {scope:?}; no open session grants that exact namespace and no flush cursor was advanced"
            ),
            SessionError::LedgerScopeSinkMismatch {
                scope,
                bound_sink,
                attempted_sink,
            } => write!(
                f,
                "ledger scope {scope:?} is already bound to ledger instance {bound_sink}; refusing instance {attempted_sink} and leaving every scope cursor unchanged"
            ),
            SessionError::ScopePermitMismatch { scope } => write!(
                f,
                "scope flush permit for {scope:?} belongs to a different governor; no sink or cursor state was changed"
            ),
            SessionError::ScopeFlushInFlight { scope } => write!(
                f,
                "ledger scope {scope:?} already has a bounded flush in flight; retry after it completes"
            ),
            SessionError::LimitExceeded {
                resource,
                limit,
                observed_at_least,
            } => write!(
                f,
                "session {resource} limit {limit} exceeded (observed at least {observed_at_least}); no partial authority mutation was committed"
            ),
            SessionError::InvalidResource {
                resource,
                value,
                requirement,
            } => write!(
                f,
                "invalid {resource} value {value}: {requirement}; session state was not mutated"
            ),
            SessionError::Submission { what } => write!(f, "submission failed: {what}"),
            SessionError::Persistence { what } => write!(f, "persistence failed: {what}"),
            SessionError::InvalidPressureLevel { level } => write!(
                f,
                "memory-pressure level {level} is outside the declared ladder 1..=3; \
                 out-of-ladder levels are refused, never clamped"
            ),
            SessionError::UngatedSession { id } => write!(
                f,
                "session {id} was opened without a cancellation gate; level-3 pressure \
                 (pause-serialize-resume) is refused — open with open_session_gated to \
                 bind the session's own gate"
            ),
            SessionError::PreRequestedGate { id } => write!(
                f,
                "session {id} supplied an already-requested cancellation gate; registration was refused so stale cancellation cannot become a new execution generation"
            ),
            SessionError::PauseRequestMismatch {
                id,
                requested_ordinal,
            } => write!(
                f,
                "pause request at ordinal {requested_ordinal} does not match session {id}'s live or replayable generation"
            ),
            SessionError::PauseAcknowledgementConflict {
                id,
                requested_ordinal,
            } => write!(
                f,
                "session {id} pause request at ordinal {requested_ordinal} was already acknowledged with different checkpoint evidence"
            ),
            SessionError::ResumeNotActivated { id, generation } => write!(
                f,
                "session {id} gate generation {generation} is ready to resume but not activated; pressure transitions remain refused"
            ),
            SessionError::ResumeAcknowledgementMismatch { id } => write!(
                f,
                "session {id} resume acknowledgement is foreign, stale, or inconsistent with the governor's current gate"
            ),
            SessionError::ResumeGateAlreadyRequested { id, generation } => write!(
                f,
                "session {id} resume gate generation {generation} was requested before activation; refusing to start work on a cancelled generation"
            ),
            SessionError::PauseAlreadyPending {
                id,
                requested_ordinal,
            } => write!(
                f,
                "session {id} already has a pause request pending at ordinal {requested_ordinal}; acknowledge it before requesting another pressure transition"
            ),
        }
    }
}

impl std::error::Error for SessionError {}

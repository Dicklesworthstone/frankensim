//! fs-session (plan §11.3): sessions, capability tokens, and the resource
//! GOVERNOR — budgets are ENFORCED, not advisory — plus the agent-proofing
//! trio: idempotency keys (a retry cannot double-spend), `estimate()` dry
//! runs (plan before you spend), and errors as GUIDANCE ("a refusal that
//! teaches is worth ten silent successes").
//! Quiescent session snapshots can also persist and automatically surface the
//! expansion program's quantitative PR-001--PR-012 register through
//! [`program_risk`].
//!
//! Layer: L6 (HELM). Threading contract: the governor's hot paths are
//! `Send + Sync` (in-memory, mutex-guarded) so enforcement and idempotency
//! survive concurrent submission storms; ledger persistence is an explicit
//! single-threaded `flush_scope_to_ledger` step because fsqlite connections are
//! `!Send` by design.

pub mod estimate;
pub mod gemm_tune;
pub mod governor;
pub mod grant;
pub mod guidance;
pub mod program_risk;
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
    Charge, DegradationEvent, DegradationStep, DurableGovernorNonce, Enforcement, FlushReport,
    Governor, MAX_CHECKPOINT_CLAIM_BYTES, MAX_DEGRADATION_EVENTS_PER_SCOPE, MAX_EVENT_PAGE_ROWS,
    MAX_FLUSH_ENCODED_BYTES, MAX_FLUSH_ROWS, MAX_IDEMPOTENCY_INPUT_BYTES,
    MAX_IDEMPOTENCY_KEYS_PER_SESSION, MAX_METER_REPORTS_PER_SESSION,
    MAX_PRESSURE_ACTIONS_PER_SESSION, MAX_RETAINED_BYTES_PER_GOVERNOR,
    MAX_RETAINED_BYTES_PER_SCOPE, MAX_RETAINED_EVIDENCE_BYTES, MAX_SESSIONS_PER_GOVERNOR,
    MAX_SESSIONS_PER_SCOPE, MeterReceipt, MeterReportId, MeterSnapshot, PauseAcknowledgement,
    PauseRequestId, PressureActionId, PressureReceipt, ResumeActivationId, ResumeActivationReceipt,
    RetainedEvidence, ScopeFlushPermit, SessionOpenId, SessionOpenReceipt, StepPhase,
    SubmissionReceipt, SubmissionRequestId, SubmitOutcome,
};
pub use grant::{
    CoreLease, CoreLeaseBook, IssuerIdentity, IssuerPolicy, MAX_ISSUER_FIELD_BYTES, NoIssuerPolicy,
    PolicyDecision, SessionGrant, mint_grant,
};
pub use guidance::Guidance;
pub use program_risk::{
    PROGRAM_RISK_REGISTER_ARTIFACT_KIND, PROGRAM_RISK_REPORT_CODEC_VERSION,
    PROGRAM_RISK_REPORT_EVENT_KIND, PROGRAM_RISK_REPORT_ID_DOMAIN,
    PROGRAM_RISK_REPORT_IDENTITY_VERSION, PROGRAM_RISK_REPORT_LOGICAL_ROWS,
    PROGRAM_RISK_REPORT_ROW_ORDER_VERSION, PROGRAM_RISK_REPORT_STATUS_TAG_VERSION,
    PROGRAM_RISK_SESSION_REPORT_ARTIFACT_KIND, ProgramRiskAlert, ProgramRiskReportDisposition,
    ProgramRiskReportId, ProgramRiskReportReceipt, ProgramRiskReportWrite,
};
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
    /// An issuer identity field is unbounded or non-canonical (aeq7).
    InvalidIssuerField {
        /// Which field.
        field: &'static str,
        /// Observed byte length.
        observed_bytes: usize,
    },
    /// The injected issuer policy refused to mint a grant (aeq7).
    GrantDenied {
        /// The policy's teaching reason.
        reason: String,
    },
    /// A grant failed issuer/fingerprint/digest verification (aeq7).
    GrantForged {
        /// The session the grant claims.
        session: u64,
    },
    /// A grant is past its admitted expiry (aeq7).
    GrantExpired {
        /// The session.
        session: u64,
        /// Admitted exclusive expiry, ledger nanoseconds.
        expiry_ns: i64,
        /// The refusing evaluation time.
        now_ns: i64,
    },
    /// A grant's revocation generation is no longer current (aeq7).
    GrantRevoked {
        /// The session.
        session: u64,
        /// Generation admitted at mint time.
        granted_generation: u64,
        /// The policy's current generation.
        current_generation: u64,
    },
    /// Execution asked for an operator the admitted grant never named
    /// (aeq7).
    UngrantedVerb {
        /// The session.
        session: u64,
        /// The refused operator.
        verb: String,
    },
    /// A core-lease acquisition would exceed the admitted concurrency
    /// (aeq7).
    CoreLeaseExceeded {
        /// The session.
        session: u64,
        /// Admitted concurrent cores.
        granted: u64,
        /// Cores already leased.
        active: u64,
        /// Cores requested.
        requested: u64,
    },
    /// A session-end snapshot was requested while caller work or a pause
    /// acknowledgement was still in flight.
    SessionNotQuiescent {
        /// Session that has not drained.
        id: u64,
        /// Exact admitted submissions still executing.
        pending_submissions: usize,
        /// Whether a pause request still awaits acknowledgement.
        pause_pending: bool,
    },
    /// This live governor is already publishing or recovering the same
    /// program-risk singleton.
    ProgramRiskReportInFlight {
        /// Domain-separated singleton authority already reserved in-process.
        id: fs_blake3::ContentHash,
    },
    /// A durable report belongs to a gate generation whose lifecycle has not
    /// yet been recovered into this governor.
    ProgramRiskReportGenerationAhead {
        /// Session whose recovered lifecycle is behind the report.
        id: u64,
        /// Generation recorded by the durable report.
        report_generation: u64,
        /// Generation currently reconstructed by lifecycle recovery.
        recovered_generation: u64,
    },
    /// A session id was registered more than once. Session identity is
    /// immutable: replacing a token would let new authority inherit old
    /// meters, pause state, and idempotency receipts.
    SessionAlreadyOpen {
        /// The duplicate id.
        id: u64,
    },
    /// An opaque mutation authority belongs to another governor, session, or
    /// immutable session-open identity.
    MutationAuthorityMismatch {
        /// Bounded authority family.
        kind: &'static str,
        /// Domain-separated identity of the rejected authority.
        id: fs_blake3::ContentHash,
    },
    /// An already-committed retry authority was reused with a different
    /// payload or execution capability.
    MutationConflict {
        /// Bounded authority family.
        kind: &'static str,
        /// Domain-separated identity whose first payload remains authoritative.
        id: fs_blake3::ContentHash,
    },
    /// An unused authority names an execution generation that is no longer
    /// current. Exact replay of an already-committed receipt is still allowed.
    StaleMutationGeneration {
        /// Bounded authority family.
        kind: &'static str,
        /// Session named by the authority.
        id: u64,
        /// Generation captured when the authority was minted.
        supplied: u64,
        /// Current live generation.
        current: u64,
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
    /// A restart-stable governor refused a fresh execution through the
    /// process-local API because no ledger was supplied for its mandatory
    /// pre-execution claim.
    DurableLedgerRequired {
        /// Mutation family.
        kind: &'static str,
        /// Restart-stable authority that must use the durable path.
        authority: fs_blake3::ContentHash,
    },
    /// Durable recovery found no terminal receipt for an authority whose
    /// caller must explicitly choose a reconciliation path.
    RecoveryRequired {
        /// Mutation family.
        kind: &'static str,
        /// Restart-stable authority.
        authority: fs_blake3::ContentHash,
    },
    /// A restarted durable governor has not reconstructed every immutable
    /// claim observed at construction, so fresh mutations remain fenced.
    DurableRecoveryIncomplete {
        /// Number of historical claims still requiring typed recovery or
        /// explicit reconciliation.
        remaining_claims: u64,
    },
    /// A durable execution claim exists without a terminal receipt. Caller
    /// work may have produced side effects, so it is never run automatically.
    IndeterminateMutation {
        /// Mutation family.
        kind: &'static str,
        /// Restart-stable authority.
        authority: fs_blake3::ContentHash,
    },
    /// A durable terminal row, receipt codec, or owned event group violated
    /// its authenticated bounded schema.
    TerminalCorrupt {
        /// Mutation family.
        kind: &'static str,
        /// Restart-stable authority.
        authority: fs_blake3::ContentHash,
        /// Bounded structural diagnosis.
        detail: String,
    },
    /// A durable receipt uses a newer fs-session codec schema.
    UnsupportedTerminalSchema {
        /// Stored schema version.
        found: u32,
        /// Highest version understood by this build.
        supported: u32,
    },
    /// Recovery was attempted against a different physical ledger.
    RecoveryLedgerMismatch {
        /// Ledger bound into the durable governor identity.
        expected: fs_ledger::LedgerInstanceId,
        /// Ledger supplied to recovery.
        attempted: fs_ledger::LedgerInstanceId,
    },
    /// A stored terminal names a different durable governor namespace.
    RecoveryGovernorMismatch {
        /// Current governor namespace.
        expected: fs_blake3::ContentHash,
        /// Namespace stored by the terminal.
        found: fs_blake3::ContentHash,
    },
    /// Causal meter recovery skipped or reordered a committed ordinal.
    RecoveryCausalGap {
        /// Session whose meter chain is discontinuous.
        session: u64,
        /// Next required commit ordinal.
        expected: u64,
        /// Stored ordinal offered to recovery.
        found: u64,
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
    /// A solver checkpoint receipt was forged, stale, cross-session, bound to
    /// another pause, or not independently verified in the supplied ledger.
    PauseCheckpointMismatch {
        /// Session named by the pending pause request.
        id: u64,
        /// Deterministic request ordinal that remained incomplete.
        requested_ordinal: i64,
        /// Stable machine-actionable refusal class.
        reason: &'static str,
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
    /// A pause acknowledgement arrived before every admitted submission in
    /// the draining generation published its terminal meter outcome.
    PauseDrainPending {
        /// Session whose gate generation cannot rotate yet.
        id: u64,
        /// Exact number of admitted submissions still executing.
        pending_submissions: usize,
    },
    /// New work was offered after the current cancellation gate entered its
    /// draining state but before a replacement generation was activated.
    SessionGateDraining {
        /// Session refusing new work.
        id: u64,
        /// Requested gate generation that is draining.
        generation: u64,
    },
}

impl fmt::Display for SessionError {
    #[allow(clippy::too_many_lines)] // Exhaustive rendering stays adjacent to the error variants.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SessionError::UnknownSession { id } => write!(f, "unknown session {id}"),
            SessionError::SessionNotQuiescent {
                id,
                pending_submissions,
                pause_pending,
            } => write!(
                f,
                "session {id} is not quiescent: {pending_submissions} submission(s) remain in flight and pause_pending={pause_pending}; drain work and acknowledge any pause before publishing a session-end snapshot"
            ),
            SessionError::ProgramRiskReportInFlight { id } => write!(
                f,
                "program-risk report {id} already has an in-process publication or recovery attempt; retry after that bounded attempt completes"
            ),
            SessionError::ProgramRiskReportGenerationAhead {
                id,
                report_generation,
                recovered_generation,
            } => write!(
                f,
                "program-risk report for session {id} belongs to generation {report_generation}, but lifecycle recovery has reconstructed only generation {recovered_generation}; recover the intervening pause/resume lifecycle before this report"
            ),
            SessionError::SessionAlreadyOpen { id } => write!(
                f,
                "session {id} is already open; capability tokens are immutable and the existing \
                 session state was left unchanged"
            ),
            SessionError::MutationAuthorityMismatch { kind, id } => write!(
                f,
                "{kind} authority {id} belongs to another governor, session, or immutable open identity; no mutation was committed"
            ),
            SessionError::MutationConflict { kind, id } => write!(
                f,
                "{kind} authority {id} was already claimed or committed with a different payload; the original authority remains authoritative"
            ),
            SessionError::StaleMutationGeneration {
                kind,
                id,
                supplied,
                current,
            } => write!(
                f,
                "{kind} authority for session {id} captured stale generation {supplied}; current generation is {current} and no mutation was committed"
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
            SessionError::InvalidIssuerField {
                field,
                observed_bytes,
            } => write!(
                f,
                "issuer {field} must be 1..={} ASCII graphic bytes, got {observed_bytes}",
                grant::MAX_ISSUER_FIELD_BYTES
            ),
            SessionError::GrantDenied { reason } => {
                write!(f, "session grant denied by the injected policy: {reason}")
            }
            SessionError::GrantForged { session } => write!(
                f,
                "grant for session {session} failed issuer/digest verification; \
                 no authority was extended"
            ),
            SessionError::GrantExpired {
                session,
                expiry_ns,
                now_ns,
            } => write!(
                f,
                "grant for session {session} expired at {expiry_ns} ns (evaluated at {now_ns} ns); \
                 request re-issuance"
            ),
            SessionError::GrantRevoked {
                session,
                granted_generation,
                current_generation,
            } => write!(
                f,
                "grant for session {session} was minted under revocation generation \
                 {granted_generation} but the policy is at {current_generation}; \
                 request re-endorsement"
            ),
            SessionError::UngrantedVerb { session, verb } => write!(
                f,
                "session {session} asked to execute operator {verb:?} outside its admitted \
                 grant; the lease was refused"
            ),
            SessionError::CoreLeaseExceeded {
                session,
                granted,
                active,
                requested,
            } => write!(
                f,
                "session {session} core lease refused: {active} active + {requested} requested \
                 exceeds the admitted {granted} concurrent cores"
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
            SessionError::DurableLedgerRequired { kind, authority } => write!(
                f,
                "durable {kind} authority {authority} requires the ledger-bound submission API before fresh caller work can run"
            ),
            SessionError::RecoveryRequired { kind, authority } => write!(
                f,
                "durable {kind} authority {authority} requires explicit typed recovery or reconciliation"
            ),
            SessionError::DurableRecoveryIncomplete { remaining_claims } => write!(
                f,
                "durable governor recovery is incomplete: {remaining_claims} historical mutation claims remain; fresh mutation is refused"
            ),
            SessionError::IndeterminateMutation { kind, authority } => write!(
                f,
                "durable {kind} authority {authority} was claimed but has no terminal receipt; caller work may have run, so automatic re-execution is refused"
            ),
            SessionError::TerminalCorrupt {
                kind,
                authority,
                detail,
            } => write!(
                f,
                "durable {kind} terminal {authority} is corrupt: {detail}; typed replay is refused"
            ),
            SessionError::UnsupportedTerminalSchema { found, supported } => write!(
                f,
                "durable terminal codec schema v{found} is newer than supported v{supported}; upgrade before recovery"
            ),
            SessionError::RecoveryLedgerMismatch {
                expected,
                attempted,
            } => write!(
                f,
                "durable recovery is bound to ledger {expected}; refusing foreign ledger {attempted}"
            ),
            SessionError::RecoveryGovernorMismatch { expected, found } => write!(
                f,
                "durable terminal belongs to governor {found}; current recovery namespace is {expected}"
            ),
            SessionError::RecoveryCausalGap {
                session,
                expected,
                found,
            } => write!(
                f,
                "durable meter recovery for session {session} expected commit ordinal {expected} but found {found}; recover the contiguous prefix first"
            ),
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
            SessionError::PauseCheckpointMismatch {
                id,
                requested_ordinal,
                reason,
            } => write!(
                f,
                "session {id} pause request at ordinal {requested_ordinal} refused checkpoint receipt ({reason}); the request remains pending"
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
            SessionError::PauseDrainPending {
                id,
                pending_submissions,
            } => write!(
                f,
                "session {id} still has {pending_submissions} admitted submission(s) draining; their terminal meter outcomes must publish before the pause generation can rotate"
            ),
            SessionError::SessionGateDraining { id, generation } => write!(
                f,
                "session {id} gate generation {generation} is draining after cancellation; exact terminal replays remain available but new work is refused"
            ),
        }
    }
}

impl std::error::Error for SessionError {}

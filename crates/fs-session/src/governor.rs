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

/// Hard-bound ratio: past 6/5 of a grant the session pauses. Float and exact
/// integer resource paths derive from this one policy definition.
const HARD_FACTOR_NUMERATOR: u32 = 6;
const HARD_FACTOR_DENOMINATOR: u32 = 5;
#[allow(clippy::cast_lossless)] // small policy integers are exactly representable as f64
const HARD_FACTOR: f64 = HARD_FACTOR_NUMERATOR as f64 / HARD_FACTOR_DENOMINATOR as f64;
const IDEMPOTENCY_KEY_DOMAIN: &str = "org.frankensim.fs-session.idempotency-key.v2";
const SUBMISSION_RECEIPT_DOMAIN: &str = "org.frankensim.fs-session.submission-receipt.v2";
const RETAINED_EVIDENCE_DOMAIN: &str = "org.frankensim.fs-session.retained-evidence.v1";
const MAX_IDEMPOTENCY_KEY_BYTES: usize = 4096;
// Fixed-width conservative framing for t, three byte lengths, and the two
// optional-field discriminants in one persisted event row.
const FLUSH_ROW_FRAMING_BYTES: usize =
    core::mem::size_of::<i64>() + 3 * core::mem::size_of::<u64>() + 2 * core::mem::size_of::<u8>();

/// Maximum sessions admitted by one governor.
pub const MAX_SESSIONS_PER_GOVERNOR: usize = 4096;
/// Maximum sessions sharing one exact ledger scope.
pub const MAX_SESSIONS_PER_SCOPE: usize = 1024;
/// Maximum distinct idempotency keys retained for one session.
pub const MAX_IDEMPOTENCY_KEYS_PER_SESSION: usize = 4096;
/// Maximum degradation events retained in memory for one scope.
pub const MAX_DEGRADATION_EVENTS_PER_SCOPE: usize = 65_536;
/// Maximum UTF-8 bytes retained from caller-controlled diagnostic evidence.
pub const MAX_RETAINED_EVIDENCE_BYTES: usize = 16 * 1024;
/// Maximum event rows emitted by one bounded flush call.
pub const MAX_FLUSH_ROWS: usize = 1024;
/// Maximum encoded event bytes emitted by one bounded flush call.
pub const MAX_FLUSH_ENCODED_BYTES: usize = 4 * 1024 * 1024;
/// Maximum degradation events returned by one page request.
pub const MAX_EVENT_PAGE_ROWS: usize = 1024;

static NEXT_GOVERNOR_ID: AtomicU64 = AtomicU64::new(1);

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

/// Unforgeable authority to flush one exact scope from one governor.
#[derive(Clone, PartialEq, Eq)]
pub struct ScopeFlushPermit {
    governor_id: u64,
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

/// Result of one bounded flush chunk.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FlushReport {
    /// Rows atomically appended by this call.
    pub appended_rows: usize,
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
    /// The step's action was applied synchronously (spill/coarsen).
    Applied,
    /// Cancellation was requested on the session's OWN gate; the solver
    /// has not yet acknowledged with a checkpoint receipt.
    Requested,
    /// The solver acknowledged: checkpoint receipt recorded.
    Complete,
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
    /// Bounded checkpoint evidence carried by completion events.
    pub checkpoint: Option<RetainedEvidence>,
}

/// Opaque content identity for one terminal idempotent submission.
///
/// The private field prevents callers from minting receipts from arbitrary
/// integers. Identity binds the owning session, exact idempotency key, terminal
/// outcome, and charge or failure diagnosis.
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
        session: SessionId,
        ledger_scope: &str,
        idem_key: &str,
        charge: Charge,
        enforcement: &Enforcement,
    ) -> bool {
        self == submission_receipt(
            session,
            ledger_scope,
            idem_key,
            &SubmissionCompletion::Done(charge, enforcement.clone()),
        )
    }

    /// Recompute and verify a failed terminal receipt.
    #[must_use]
    pub fn matches_failure(
        self,
        session: SessionId,
        ledger_scope: &str,
        idem_key: &str,
        evidence: &RetainedEvidence,
    ) -> bool {
        self == submission_receipt(
            session,
            ledger_scope,
            idem_key,
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
        /// The charge recorded.
        charge: Charge,
        /// Enforcement decision produced by committing that charge.
        enforcement: Enforcement,
        /// Content-derived terminal receipt.
        receipt: SubmissionReceipt,
    },
    /// The key had already executed (or raced and lost): same receipt,
    /// NO additional charge.
    Duplicate {
        /// The original execution's receipt.
        receipt: SubmissionReceipt,
        /// The original execution's enforcement decision.
        enforcement: Enforcement,
    },
    /// The one attempted execution failed before a charge could be committed.
    /// The key remains terminal: all duplicates receive this same receipt and
    /// diagnosis, and an explicit retry requires a new key.
    Failed {
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

fn same_meter_snapshot(left: &SessionMeters, right: &SessionMeters) -> bool {
    left.core_s.to_bits() == right.core_s.to_bits()
        && left.mem_peak_bytes == right.mem_peak_bytes
        && left.wall_s.to_bits() == right.wall_s.to_bits()
        && left.throttled == right.throttled
        && left.paused == right.paused
}

#[derive(Debug)]
enum IdemState {
    Pending,
    Done {
        ordinal: u64,
        receipt: SubmissionReceipt,
        charge: Charge,
        enforcement: Enforcement,
    },
    Failed {
        ordinal: u64,
        receipt: SubmissionReceipt,
        evidence: RetainedEvidence,
    },
}

enum SubmissionCompletion {
    Done(Charge, Enforcement),
    Failed(RetainedEvidence),
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
    session: SessionId,
    ledger_scope: &str,
    idem_key: &str,
    completion: &SubmissionCompletion,
) -> SubmissionReceipt {
    let mut payload = Vec::new();
    payload.extend_from_slice(&session.0.to_le_bytes());
    push_framed(&mut payload, ledger_scope.as_bytes());
    push_framed(&mut payload, idem_key.as_bytes());
    match completion {
        SubmissionCompletion::Done(charge, enforcement) => {
            payload.push(0);
            payload.extend_from_slice(&charge.core_s.to_bits().to_le_bytes());
            payload.extend_from_slice(&charge.mem_peak_bytes.to_le_bytes());
            payload.extend_from_slice(&charge.wall_s.to_bits().to_le_bytes());
            push_enforcement_identity(&mut payload, enforcement);
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

struct PreparedFlush {
    reservation_id: u64,
    generation: i64,
    revision: u64,
    next_flush_lane: u8,
    buffered: Vec<BufferedLedgerEvent>,
    encoded_bytes: usize,
    meter_marks: Vec<(u64, SessionMeters)>,
    idempotency_marks: Vec<((u64, String), (u64, SubmissionReceipt))>,
    event_target: usize,
}

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
    dirty_meters: BTreeSet<u64>,
    dirty_idempotency: BTreeSet<(u64, String)>,
    events: Vec<DegradationEvent>,
    flushed_events: usize,
    sink: Option<fs_ledger::LedgerInstanceId>,
    flush_generation: i64,
    in_flight: Option<u64>,
    revision: u64,
    next_flush_lane: u8,
}

#[derive(Default)]
struct Inner {
    tokens: BTreeMap<u64, CapabilityToken>,
    /// Session-OWNED cancellation gates, bound at open (gp3.13): the
    /// only route to a pause request, so a foreign gate is
    /// unrepresentable at the pressure API.
    gates: BTreeMap<u64, Arc<CancelGate>>,
    /// Pause requests awaiting a checkpoint acknowledgement, keyed by
    /// session → ordinal of the Requested event.
    pending_pause: BTreeMap<u64, i64>,
    meters: BTreeMap<u64, SessionMeters>,
    idempotency: BTreeMap<(u64, String), IdemState>,
    idempotency_keys: BTreeMap<u64, BTreeSet<String>>,
    scopes: BTreeMap<String, ScopeState>,
    next_submission_ordinal: u64,
    next_ordinal: i64,
    next_flush_reservation: u64,
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
    id: u64,
    inner: Mutex<Inner>,
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
            id: NEXT_GOVERNOR_ID.fetch_add(1, Ordering::Relaxed),
            inner: Mutex::new(Inner::default()),
        }
    }

    fn register_session(
        &self,
        token: CapabilityToken,
        gate: Option<Arc<CancelGate>>,
    ) -> Result<ScopeFlushPermit, SessionError> {
        CapabilityToken::validate_ledger_scope(&token.ledger_scope)?;
        validate_resource("core-seconds grant", token.core_s)?;
        validate_resource("wall-seconds grant", token.wall_s)?;
        let session = token.session.0;
        let ledger_scope = token.ledger_scope.clone();
        let mut g = self.inner.lock().expect("governor lock");
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
        let next_revision = g
            .scopes
            .get(&ledger_scope)
            .map_or(1, |scope| scope.revision.saturating_add(1));
        g.meters.insert(session, SessionMeters::default());
        g.idempotency_keys.insert(session, BTreeSet::new());
        g.tokens.insert(session, token);
        if let Some(gate) = gate {
            g.gates.insert(session, gate);
        }
        let scope = g.scopes.entry(ledger_scope.clone()).or_default();
        scope.sessions.insert(session);
        scope.dirty_meters.insert(session);
        scope.revision = next_revision;
        Ok(ScopeFlushPermit {
            governor_id: self.id,
            ledger_scope,
        })
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
    pub fn open_session(&self, token: CapabilityToken) -> Result<ScopeFlushPermit, SessionError> {
        self.register_session(token, None)
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
        token: CapabilityToken,
        gate: Arc<CancelGate>,
    ) -> Result<ScopeFlushPermit, SessionError> {
        self.register_session(token, Some(gate))
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
            .cloned()
            .ok_or(SessionError::UnknownSession { id: session.0 })
    }

    /// Meter a consumption delta and enforce the token bounds:
    /// at the grant → `Throttled`; past `HARD_FACTOR ×` grant → `Paused`.
    /// Structured outcomes only — the governor NEVER silently kills.
    ///
    /// # Errors
    /// [`SessionError::UnknownSession`], [`SessionError::InvalidResource`], or
    /// [`SessionError::LimitExceeded`] if a throttle/pause counter is exhausted.
    pub fn charge(&self, session: SessionId, delta: Charge) -> Result<Enforcement, SessionError> {
        validate_resource("core-seconds charge", delta.core_s)?;
        validate_resource("wall-seconds charge", delta.wall_s)?;
        let mut g = self.inner.lock().expect("governor lock");
        let token = g
            .tokens
            .get(&session.0)
            .cloned()
            .ok_or(SessionError::UnknownSession { id: session.0 })?;
        let mut next = g.meters.get(&session.0).cloned().unwrap_or_default();
        let next_core_s = next.core_s + delta.core_s;
        let next_wall_s = next.wall_s + delta.wall_s;
        validate_resource("accumulated core-seconds", next_core_s)?;
        validate_resource("accumulated wall-seconds", next_wall_s)?;
        next.core_s = next_core_s;
        next.mem_peak_bytes = next.mem_peak_bytes.max(delta.mem_peak_bytes);
        next.wall_s = next_wall_s;
        // Memory is an exact byte budget. Converting u64 values to f64 before
        // admission collapses adjacent values above 2^53 and can throttle a
        // session that is still below its grant. Compare the 6/5 hard boundary
        // exactly in u128; f64 remains only the legacy diagnostic projection.
        let memory_past_hard = u128::from(next.mem_peak_bytes)
            * u128::from(HARD_FACTOR_DENOMINATOR)
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
                    "checkpoint required before continuing; resume with a larger {resource} \
                     grant or a coarsened study — the caller must arrange and ledger the \
                     checkpoint explicitly"
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
                next.throttled =
                    next.throttled
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
        bump_scope_revision(&mut g, &token.ledger_scope);
        g.meters.insert(session.0, next);
        g.scopes
            .get_mut(&token.ledger_scope)
            .expect("registered session scope")
            .dirty_meters
            .insert(session.0);
        Ok(enforcement)
    }

    /// Idempotency-keyed exactly-once execution: the first caller runs `work`
    /// and is charged. A concurrent caller observes [`SubmitOutcome::InFlight`]
    /// without blocking; after terminal publication, repeat callers receive the
    /// same receipt and no additional charge.
    ///
    /// # Errors
    /// [`SessionError::UnknownSession`] for an unknown owner,
    /// [`SessionError::Submission`] for a blank/oversized key, or
    /// [`SessionError::LimitExceeded`] for retained-key or logical-ordinal
    /// exhaustion.
    ///
    /// A panic in `work` is contained and committed as a terminal
    /// [`SubmitOutcome::Failed`] receipt. The same key never reruns implicitly:
    /// duplicates receive that same failure receipt and callers must choose a
    /// new idempotency key for an explicit retry.
    #[allow(clippy::too_many_lines)] // Keep the Pending-to-terminal transaction contiguous.
    pub fn submit_once(
        &self,
        session: SessionId,
        idem_key: &str,
        work: impl FnOnce() -> Charge,
    ) -> Result<SubmitOutcome, SessionError> {
        if idem_key.trim().is_empty() || idem_key.len() > MAX_IDEMPOTENCY_KEY_BYTES {
            return Err(SessionError::Submission {
                what: format!(
                    "idempotency key must be non-blank and at most {MAX_IDEMPOTENCY_KEY_BYTES} bytes"
                ),
            });
        }
        let idempotency_scope = (session.0, idem_key.to_string());
        let (ordinal, ledger_scope) = {
            let mut g = self.inner.lock().expect("governor lock");
            let ledger_scope = g
                .tokens
                .get(&session.0)
                .map(|token| token.ledger_scope.clone())
                .ok_or(SessionError::UnknownSession { id: session.0 })?;
            match g.idempotency.get(&idempotency_scope) {
                Some(IdemState::Done {
                    receipt,
                    enforcement,
                    ..
                }) => {
                    return Ok(SubmitOutcome::Duplicate {
                        receipt: *receipt,
                        enforcement: enforcement.clone(),
                    });
                }
                Some(IdemState::Failed {
                    receipt, evidence, ..
                }) => {
                    return Ok(SubmitOutcome::Failed {
                        receipt: *receipt,
                        evidence: evidence.clone(),
                    });
                }
                Some(IdemState::Pending) => return Ok(SubmitOutcome::InFlight),
                None => {}
            }
            let key_count = g.idempotency_keys.get(&session.0).map_or(0, BTreeSet::len);
            if key_count >= MAX_IDEMPOTENCY_KEYS_PER_SESSION {
                return Err(SessionError::LimitExceeded {
                    resource: "idempotency_keys_per_session",
                    limit: MAX_IDEMPOTENCY_KEYS_PER_SESSION,
                    observed_at_least: key_count.saturating_add(1),
                });
            }
            let ordinal =
                g.next_submission_ordinal
                    .checked_add(1)
                    .ok_or(SessionError::LimitExceeded {
                        resource: "submission_ordinal",
                        limit: i64::MAX as usize,
                        observed_at_least: usize::MAX,
                    })?;
            if ordinal > i64::MAX as u64 {
                return Err(SessionError::LimitExceeded {
                    resource: "submission_ordinal",
                    limit: i64::MAX as usize,
                    observed_at_least: i64::MAX as usize,
                });
            }
            g.next_submission_ordinal = ordinal;
            g.idempotency
                .insert(idempotency_scope.clone(), IdemState::Pending);
            g.idempotency_keys
                .get_mut(&session.0)
                .expect("registered session key index")
                .insert(idem_key.to_string());
            (ordinal, ledger_scope)
        };
        // Execute OUTSIDE the lock (work may be long). Catching here is
        // load-bearing: every Pending key must reach a terminal state even when
        // caller-authored work unwinds.
        let completion = match std::panic::catch_unwind(std::panic::AssertUnwindSafe(work)) {
            Ok(charge) => match self.charge(session, charge) {
                Ok(enforcement) => SubmissionCompletion::Done(charge, enforcement),
                Err(error) => {
                    SubmissionCompletion::Failed(RetainedEvidence::capture(&error.to_string()))
                }
            },
            Err(payload) => SubmissionCompletion::Failed(panic_evidence(payload.as_ref())),
        };
        let receipt = submission_receipt(session, &ledger_scope, idem_key, &completion);
        let outcome;
        {
            let mut g = self.inner.lock().expect("governor lock");
            let dirty_idempotency = idempotency_scope.clone();
            match completion {
                SubmissionCompletion::Done(charge, enforcement) => {
                    g.idempotency.insert(
                        idempotency_scope,
                        IdemState::Done {
                            ordinal,
                            receipt,
                            charge,
                            enforcement: enforcement.clone(),
                        },
                    );
                    outcome = SubmitOutcome::Executed {
                        charge,
                        enforcement,
                        receipt,
                    };
                }
                SubmissionCompletion::Failed(evidence) => {
                    g.idempotency.insert(
                        idempotency_scope,
                        IdemState::Failed {
                            ordinal,
                            receipt,
                            evidence: evidence.clone(),
                        },
                    );
                    outcome = SubmitOutcome::Failed { receipt, evidence };
                }
            }
            bump_scope_revision(&mut g, &ledger_scope);
            g.scopes
                .get_mut(&ledger_scope)
                .expect("registered session scope")
                .dirty_idempotency
                .insert(dirty_idempotency);
        }
        Ok(outcome)
    }

    /// The canonical idempotency key: length-framed agent key and program text
    /// under a domain-separated BLAKE3 identity.
    #[must_use]
    pub fn idempotency_key(agent_key: &str, program_text: &str) -> String {
        let mut payload = Vec::new();
        push_framed(&mut payload, agent_key.as_bytes());
        push_framed(&mut payload, program_text.as_bytes());
        format!(
            "fs-session-idem-v2:{}",
            fs_blake3::hash_domain(IDEMPOTENCY_KEY_DOMAIN, &payload)
        )
    }

    /// Apply memory pressure at `level` (1..=3 ONLY): ladder steps
    /// `1..=level` fire IN THE DECLARED ORDER, each recorded with
    /// attribution. The `PauseSerializeResume` step requests
    /// cancellation on the session's OWN gate, resolved by `SessionId`
    /// from the binding made at [`Governor::open_session_gated`] — no
    /// gate crosses this API, so pausing a different session's work is
    /// unrepresentable (bead gp3.13). The request event is phase
    /// `Requested`; it becomes `Complete` only through
    /// [`Governor::acknowledge_pause`] with a checkpoint receipt.
    ///
    /// # Errors
    /// - [`SessionError::InvalidPressureLevel`] for levels 0 and > 3.
    /// - [`SessionError::UnknownSession`].
    /// - [`SessionError::UngatedSession`] when level 3 targets a
    ///   session opened without a cancellation gate. Validation is
    ///   ATOMIC: no step fires and nothing is ledgered.
    /// - [`SessionError::PauseAlreadyPending`] when level 3 would overwrite an
    ///   unacknowledged request.
    /// - [`SessionError::LimitExceeded`] for event or ordinal exhaustion.
    #[allow(clippy::too_many_lines)] // The ordered preflight and ladder commit are one state machine.
    pub fn apply_memory_pressure(
        &self,
        session: SessionId,
        level: u8,
    ) -> Result<Vec<DegradationEvent>, SessionError> {
        if !(1..=3).contains(&level) {
            return Err(SessionError::InvalidPressureLevel { level });
        }
        let mut g = self.inner.lock().expect("governor lock");
        let ledger_scope = g
            .tokens
            .get(&session.0)
            .map(|token| token.ledger_scope.clone())
            .ok_or(SessionError::UnknownSession { id: session.0 })?;
        if usize::from(level) >= LADDER.len()
            && let Some(requested_ordinal) = g.pending_pause.get(&session.0)
        {
            return Err(SessionError::PauseAlreadyPending {
                id: session.0,
                requested_ordinal: *requested_ordinal,
            });
        }
        let scope_event_count = g
            .scopes
            .get(&ledger_scope)
            .expect("registered session scope")
            .events
            .len();
        let requested_event_count = scope_event_count.saturating_add(usize::from(level));
        if requested_event_count > MAX_DEGRADATION_EVENTS_PER_SCOPE {
            return Err(SessionError::LimitExceeded {
                resource: "degradation_events_per_scope",
                limit: MAX_DEGRADATION_EVENTS_PER_SCOPE,
                observed_at_least: requested_event_count,
            });
        }
        let final_ordinal =
            g.next_ordinal
                .checked_add(i64::from(level))
                .ok_or(SessionError::LimitExceeded {
                    resource: "degradation_ordinal",
                    limit: i64::MAX as usize,
                    observed_at_least: usize::MAX,
                })?;
        // Resolve the session's own gate BEFORE any step fires: a
        // refused level-3 request must not half-apply the ladder.
        let gate = if usize::from(level) >= LADDER.len() {
            Some(
                g.gates
                    .get(&session.0)
                    .cloned()
                    .ok_or(SessionError::UngatedSession { id: session.0 })?,
            )
        } else {
            None
        };
        if let Some(gate) = &gate {
            gate.request();
        }
        let first_ordinal = g.next_ordinal + 1;
        let mut fired = Vec::with_capacity(usize::from(level));
        for (i, step) in LADDER.iter().enumerate() {
            if i >= usize::from(level) {
                break;
            }
            let (phase, attribution) = match step {
                DegradationStep::SpillColdArenas => (
                    StepPhase::Applied,
                    "spilled coldest arenas (least-recently-touched first)".to_string(),
                ),
                DegradationStep::CoarsenAdaptively => (
                    StepPhase::Applied,
                    "coarsened adaptive resolutions outside protected bands".to_string(),
                ),
                DegradationStep::PauseSerializeResume => (
                    StepPhase::Requested,
                    "requested pause on the session-owned gate: solver checkpoints \
                         at the next tile boundary (SolverState snapshot to the ledger); \
                         complete only on acknowledge_pause with a checkpoint receipt"
                        .to_string(),
                ),
            };
            let event = DegradationEvent {
                session,
                step: *step,
                pressure_level: level,
                phase,
                attribution,
                ordinal: first_ordinal
                    + i64::try_from(i).expect("the fixed degradation ladder length fits i64"),
                requested_ordinal: None,
                checkpoint: None,
            };
            fired.push(event.clone());
        }
        g.next_ordinal = final_ordinal;
        if let Some(requested) = fired
            .iter()
            .find(|event| event.phase == StepPhase::Requested)
        {
            g.pending_pause.insert(session.0, requested.ordinal);
        }
        g.scopes
            .get_mut(&ledger_scope)
            .expect("registered session scope")
            .events
            .extend(fired.iter().cloned());
        bump_scope_revision(&mut g, &ledger_scope);
        Ok(fired)
    }

    /// Acknowledge a pending pause with the solver's checkpoint receipt
    /// (bead gp3.13): the ONLY route to a `Complete` pause event. A
    /// pause that was never requested, or a blank receipt, is refused —
    /// a missing acknowledgement can never be ledgered as complete.
    ///
    /// # Errors
    /// - [`SessionError::UnknownSession`].
    /// - [`SessionError::Submission`] for a blank checkpoint receipt
    ///   (refused BEFORE the pending request is consumed).
    /// - [`SessionError::NoPendingPause`] when no pause request is
    ///   outstanding for the session.
    /// - [`SessionError::LimitExceeded`] for event or ordinal exhaustion.
    pub fn acknowledge_pause(
        &self,
        session: SessionId,
        checkpoint_receipt: &str,
    ) -> Result<DegradationEvent, SessionError> {
        if checkpoint_receipt.trim().is_empty() {
            return Err(SessionError::Submission {
                what: "pause acknowledgement requires a non-empty checkpoint receipt".to_string(),
            });
        }
        let evidence = RetainedEvidence::capture(checkpoint_receipt);
        let mut g = self.inner.lock().expect("governor lock");
        let ledger_scope = g
            .tokens
            .get(&session.0)
            .map(|token| token.ledger_scope.clone())
            .ok_or(SessionError::UnknownSession { id: session.0 })?;
        let requested_ordinal = *g
            .pending_pause
            .get(&session.0)
            .ok_or(SessionError::NoPendingPause { id: session.0 })?;
        let event_count = g
            .scopes
            .get(&ledger_scope)
            .expect("registered session scope")
            .events
            .len();
        if event_count >= MAX_DEGRADATION_EVENTS_PER_SCOPE {
            return Err(SessionError::LimitExceeded {
                resource: "degradation_events_per_scope",
                limit: MAX_DEGRADATION_EVENTS_PER_SCOPE,
                observed_at_least: event_count.saturating_add(1),
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
                "pause complete: checkpoint evidence {:?} ({} bytes, digest {}) acknowledges \
                 the request at ordinal {requested_ordinal}",
                evidence.preview, evidence.byte_len, evidence.digest
            ),
            ordinal,
            requested_ordinal: Some(requested_ordinal),
            checkpoint: Some(evidence),
        };
        g.pending_pause.remove(&session.0);
        g.next_ordinal = ordinal;
        g.scopes
            .get_mut(&ledger_scope)
            .expect("registered session scope")
            .events
            .push(event.clone());
        bump_scope_revision(&mut g, &ledger_scope);
        Ok(event)
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
        Ok(events[start..end].to_vec())
    }

    /// Persist at most one bounded chunk for the permit's exact scope.
    /// Preparation and cursor commit hold the governor mutex; atomic database
    /// I/O does not. Call again while [`FlushReport::remaining_dirty`] is true.
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
        let sink_identity = ledger.instance_id();
        let prepared = {
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
            let start_flush_lane = scope.next_flush_lane;
            let next_flush_lane = (start_flush_lane + 1) % 3;
            let mut event_target = scope.flushed_events;
            let event_count = scope.events.len();
            if event_target > event_count {
                return Err(SessionError::Persistence {
                    what: format!(
                        "scope event cursor {event_target} exceeds event count {event_count}"
                    ),
                });
            }

            let mut buffered = Vec::with_capacity(MAX_FLUSH_ROWS.min(64));
            let mut encoded_bytes = 0usize;
            let mut meter_marks = Vec::new();
            let mut idempotency_marks = Vec::new();

            // Rotate the first lane after every successful non-empty chunk.
            // This bounds preparation by dirty rows rather than retained state
            // and prevents continuously dirty meters from starving terminal
            // receipts or degradation events.
            'lanes: for lane_offset in 0..3 {
                let remaining_rows = MAX_FLUSH_ROWS - buffered.len();
                if remaining_rows == 0 {
                    break;
                }
                match (start_flush_lane + lane_offset) % 3 {
                    0 => {
                        let (dirty, has_more) = {
                            let scope = g.scopes.get(&ledger_scope).expect("scope checked above");
                            let dirty: Vec<u64> = scope
                                .dirty_meters
                                .iter()
                                .take(remaining_rows)
                                .copied()
                                .collect();
                            let has_more = scope.dirty_meters.len() > dirty.len();
                            (dirty, has_more)
                        };
                        for id in dirty {
                            let meters =
                                g.meters.get(&id).ok_or_else(|| SessionError::Persistence {
                                    what: format!(
                                        "scope dirty-meter index references missing session {id}"
                                    ),
                                })?;
                            let event = BufferedLedgerEvent {
                                session: id.to_be_bytes(),
                                t: generation,
                                kind: "session.consumption",
                                payload: scoped_payload(
                                    "fs-session-consumption-v3",
                                    &ledger_scope,
                                    &format!(
                                        "\"flush_generation\":{generation},\"core_s\":{},\"mem_peak\":{},\"wall_s\":{},\"throttled\":{},\"paused\":{}",
                                        meters.core_s,
                                        meters.mem_peak_bytes,
                                        meters.wall_s,
                                        meters.throttled,
                                        meters.paused,
                                    ),
                                ),
                            };
                            if !push_bounded_event(&mut buffered, &mut encoded_bytes, event)? {
                                break 'lanes;
                            }
                            meter_marks.push((id, meters.clone()));
                        }
                        if has_more {
                            break 'lanes;
                        }
                    }
                    1 => {
                        let (dirty, has_more) = {
                            let scope = g.scopes.get(&ledger_scope).expect("scope checked above");
                            let dirty: Vec<(u64, String)> = scope
                                .dirty_idempotency
                                .iter()
                                .take(remaining_rows)
                                .cloned()
                                .collect();
                            let has_more = scope.dirty_idempotency.len() > dirty.len();
                            (dirty, has_more)
                        };
                        for idempotency_scope in dirty {
                            let (session, key) = &idempotency_scope;
                            let state = g.idempotency.get(&idempotency_scope).ok_or_else(|| {
                                SessionError::Persistence {
                                    what: format!(
                                        "scope dirty-idempotency index references missing session {session} key"
                                    ),
                                }
                            })?;
                            let (ordinal, receipt, kind, body) = match state {
                                IdemState::Pending => {
                                    return Err(SessionError::Persistence {
                                        what: format!(
                                            "scope dirty-idempotency index references pending session {session} key"
                                        ),
                                    });
                                }
                                IdemState::Done {
                                    ordinal,
                                    receipt,
                                    charge,
                                    enforcement,
                                } => (
                                    *ordinal,
                                    *receipt,
                                    "session.idempotent-execution",
                                    format!(
                                        "\"session\":{session},\"key\":\"{}\",\"receipt\":\"{receipt}\",\"core_s_bits\":\"{:016x}\",\"mem_peak_bytes\":{},\"wall_s_bits\":\"{:016x}\",\"enforcement\":{}",
                                        json_escape(key),
                                        charge.core_s.to_bits(),
                                        charge.mem_peak_bytes,
                                        charge.wall_s.to_bits(),
                                        enforcement_json(enforcement),
                                    ),
                                ),
                                IdemState::Failed {
                                    ordinal,
                                    receipt,
                                    evidence,
                                } => (
                                    *ordinal,
                                    *receipt,
                                    "session.idempotent-failure",
                                    format!(
                                        "\"session\":{session},\"key\":\"{}\",\"receipt\":\"{receipt}\",\"error_evidence\":{}",
                                        json_escape(key),
                                        evidence_json(evidence),
                                    ),
                                ),
                            };
                            let terminal_generation = (ordinal, receipt);
                            let event = BufferedLedgerEvent {
                                session: session.to_be_bytes(),
                                t: i64::try_from(ordinal).map_err(|_| {
                                    SessionError::LimitExceeded {
                                        resource: "submission_ordinal",
                                        limit: i64::MAX as usize,
                                        observed_at_least: usize::MAX,
                                    }
                                })?,
                                kind,
                                payload: scoped_payload(
                                    "fs-session-idempotency-v4",
                                    &ledger_scope,
                                    &body,
                                ),
                            };
                            if !push_bounded_event(&mut buffered, &mut encoded_bytes, event)? {
                                break 'lanes;
                            }
                            idempotency_marks.push((idempotency_scope, terminal_generation));
                        }
                        if has_more {
                            break 'lanes;
                        }
                    }
                    2 => {
                        let scope = g.scopes.get(&ledger_scope).expect("scope checked above");
                        let end = event_target.saturating_add(remaining_rows).min(event_count);
                        for event in &scope.events[event_target..end] {
                            let requested = event
                                .requested_ordinal
                                .map_or_else(|| "null".to_string(), |ordinal| ordinal.to_string());
                            let checkpoint = event
                                .checkpoint
                                .as_ref()
                                .map_or_else(|| "null".to_string(), evidence_json);
                            let row = BufferedLedgerEvent {
                                session: event.session.0.to_be_bytes(),
                                t: event.ordinal,
                                kind: "session.degradation",
                                payload: scoped_payload(
                                    "fs-session-degradation-v3",
                                    &ledger_scope,
                                    &format!(
                                        "\"step\":\"{:?}\",\"level\":{},\"phase\":\"{:?}\",\"attribution\":\"{}\",\"requested_ordinal\":{requested},\"checkpoint\":{checkpoint}",
                                        event.step,
                                        event.pressure_level,
                                        event.phase,
                                        json_escape(&event.attribution),
                                    ),
                                ),
                            };
                            if !push_bounded_event(&mut buffered, &mut encoded_bytes, row)? {
                                break 'lanes;
                            }
                            event_target =
                                event_target
                                    .checked_add(1)
                                    .ok_or(SessionError::LimitExceeded {
                                        resource: "scope_event_cursor",
                                        limit: usize::MAX,
                                        observed_at_least: usize::MAX,
                                    })?;
                        }
                        if event_target < event_count {
                            break 'lanes;
                        }
                    }
                    _ => unreachable!("flush lane modulo three"),
                }
            }

            if buffered.is_empty() {
                return Ok(FlushReport {
                    appended_rows: 0,
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
            PreparedFlush {
                reservation_id,
                generation,
                revision,
                next_flush_lane,
                buffered,
                encoded_bytes,
                meter_marks,
                idempotency_marks,
                event_target,
            }
        };

        let rows: Vec<_> = prepared
            .buffered
            .iter()
            .map(BufferedLedgerEvent::as_row)
            .collect();
        if let Err(error) = ledger.append_events(&rows) {
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

        let appended_rows = prepared.buffered.len();
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
            scope.flushed_events = prepared.event_target;
            scope.next_flush_lane = prepared.next_flush_lane;
        }
        for (session, meters) in prepared.meter_marks {
            let still_current = g
                .meters
                .get(&session)
                .is_some_and(|current| same_meter_snapshot(current, &meters));
            if still_current {
                g.scopes
                    .get_mut(&ledger_scope)
                    .expect("committed scope remains registered")
                    .dirty_meters
                    .remove(&session);
            }
        }
        for (idempotency_scope, generation) in prepared.idempotency_marks {
            let still_current = match g.idempotency.get(&idempotency_scope) {
                Some(
                    IdemState::Done {
                        ordinal, receipt, ..
                    }
                    | IdemState::Failed {
                        ordinal, receipt, ..
                    },
                ) => (*ordinal, *receipt) == generation,
                Some(IdemState::Pending) | None => false,
            };
            if still_current {
                g.scopes
                    .get_mut(&ledger_scope)
                    .expect("committed scope remains registered")
                    .dirty_idempotency
                    .remove(&idempotency_scope);
            }
        }
        let scope = g
            .scopes
            .get(&ledger_scope)
            .expect("committed scope remains registered");
        let remaining_dirty = scope.revision != prepared.revision
            || !scope.dirty_meters.is_empty()
            || !scope.dirty_idempotency.is_empty()
            || scope.flushed_events < scope.events.len();
        Ok(FlushReport {
            appended_rows,
            encoded_bytes: prepared.encoded_bytes,
            remaining_dirty,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn retained_evidence_bounds_preview_but_receipts_bind_the_full_tail() {
        let shared_prefix = "x".repeat(MAX_RETAINED_EVIDENCE_BYTES);
        let evidence_a = RetainedEvidence::capture(&format!("{shared_prefix}A"));
        let evidence_b = RetainedEvidence::capture(&format!("{shared_prefix}B"));
        assert_eq!(evidence_a.preview(), shared_prefix);
        assert_eq!(evidence_a.preview(), evidence_b.preview());
        assert_eq!(evidence_a.byte_len(), MAX_RETAINED_EVIDENCE_BYTES + 1);
        assert_ne!(evidence_a.digest(), evidence_b.digest());

        let receipt_a = submission_receipt(
            SessionId(1),
            "scope",
            "key",
            &SubmissionCompletion::Failed(evidence_a),
        );
        let receipt_b = submission_receipt(
            SessionId(1),
            "scope",
            "key",
            &SubmissionCompletion::Failed(evidence_b),
        );
        assert_ne!(
            receipt_a, receipt_b,
            "equal retained previews must not collapse distinct full evidence"
        );
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
            phase: StepPhase::Applied,
            attribution: String::new(),
            ordinal: 0,
            requested_ordinal: None,
            checkpoint: None,
        };
        {
            let mut inner = governor.inner.lock().expect("governor lock");
            inner
                .scopes
                .get_mut("bounded")
                .expect("fixture scope")
                .events = vec![fixture; MAX_DEGRADATION_EVENTS_PER_SCOPE - 1];
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
            .submit_once(SessionId(0), "terminal", Charge::default)
            .expect("terminal fixture");
        governor
            .apply_memory_pressure(SessionId(0), 1)
            .expect("event fixture");
        let ledger = fs_ledger::Ledger::open(":memory:").expect("fixture ledger");
        let permit = permit.expect("at least one fixture session");
        let first = governor
            .flush_scope_to_ledger(&permit, &ledger)
            .expect("meter-first chunk");
        assert_eq!(first.appended_rows, MAX_FLUSH_ROWS);
        assert!(first.remaining_dirty);

        // Keep every meter dirty. The next chunk must rotate to terminal and
        // event lanes before consuming its remaining row budget on meters.
        for session in 0..MAX_SESSIONS_PER_SCOPE {
            governor
                .charge(
                    SessionId(u64::try_from(session).expect("fixture id fits")),
                    Charge::default(),
                )
                .expect("re-dirty meter");
        }
        let second = governor
            .flush_scope_to_ledger(&permit, &ledger)
            .expect("rotated chunk");
        assert_eq!(second.appended_rows, MAX_FLUSH_ROWS);
        assert!(second.remaining_dirty);
        let inner = governor.inner.lock().expect("governor lock");
        let scope = &inner.scopes["fair"];
        assert!(scope.dirty_idempotency.is_empty());
        assert_eq!(scope.flushed_events, scope.events.len());
        assert_eq!(scope.dirty_meters.len(), 2);
    }
}

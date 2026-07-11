//! The resource GOVERNOR: continuous metering against capability tokens
//! (throttle at the grant, pause past the hard bound — NEVER a silent
//! kill), idempotency-keyed exactly-once submission, and the DECLARED
//! degradation ladder under memory pressure (spill coldest arenas →
//! coarsen adaptively → pause-serialize-resume), every event recorded
//! with attribution and flushable to the Design Ledger.

use crate::token::{CapabilityToken, SessionId};
use crate::{Guidance, SessionError};
use fs_exec::CancelGate;
use std::collections::BTreeMap;
use std::sync::{Arc, Condvar, Mutex};

/// Hard-bound multiplier: past `HARD_FACTOR × grant` the session pauses.
const HARD_FACTOR: f64 = 1.2;
const IDEMPOTENCY_KEY_DOMAIN: &str = "org.frankensim.fs-session.idempotency-key.v2";
const SUBMISSION_RECEIPT_DOMAIN: &str = "org.frankensim.fs-session.submission-receipt.v1";
const MAX_IDEMPOTENCY_KEY_BYTES: usize = 4096;

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
        idem_key: &str,
        charge: Charge,
        enforcement: &Enforcement,
    ) -> bool {
        self == submission_receipt(
            session,
            idem_key,
            &SubmissionCompletion::Done(charge, enforcement.clone()),
        )
    }

    /// Recompute and verify a failed terminal receipt.
    #[must_use]
    pub fn matches_failure(self, session: SessionId, idem_key: &str, what: &str) -> bool {
        self == submission_receipt(
            session,
            idem_key,
            &SubmissionCompletion::Failed(what.to_string()),
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
        /// Panic payload or structured validation diagnosis.
        what: String,
    },
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

fn ledger_sink_identity(ledger: &fs_ledger::Ledger) -> String {
    if ledger.path() == ":memory:" {
        format!(":memory:@{ledger:p}")
    } else {
        ledger.path().to_string()
    }
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
        what: String,
    },
}

enum SubmissionCompletion {
    Done(Charge, Enforcement),
    Failed(String),
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

fn panic_message(payload: &(dyn std::any::Any + Send)) -> String {
    payload
        .downcast_ref::<&str>()
        .map(ToString::to_string)
        .or_else(|| payload.downcast_ref::<String>().cloned())
        .unwrap_or_else(|| "submission work panicked with a non-string payload".to_string())
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
    idem_key: &str,
    completion: &SubmissionCompletion,
) -> SubmissionReceipt {
    let mut payload = Vec::new();
    payload.extend_from_slice(&session.0.to_le_bytes());
    push_framed(&mut payload, idem_key.as_bytes());
    match completion {
        SubmissionCompletion::Done(charge, enforcement) => {
            payload.push(0);
            payload.extend_from_slice(&charge.core_s.to_bits().to_le_bytes());
            payload.extend_from_slice(&charge.mem_peak_bytes.to_le_bytes());
            payload.extend_from_slice(&charge.wall_s.to_bits().to_le_bytes());
            push_enforcement_identity(&mut payload, enforcement);
        }
        SubmissionCompletion::Failed(what) => {
            payload.push(1);
            push_framed(&mut payload, what.as_bytes());
        }
    }
    SubmissionReceipt(fs_blake3::hash_domain(SUBMISSION_RECEIPT_DOMAIN, &payload))
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
    events: Vec<DegradationEvent>,
    /// Last meter snapshot durably appended to the owning ledger.
    flushed_meters: BTreeMap<u64, SessionMeters>,
    /// Terminal submission generation and receipt durably appended for each
    /// session/key. Binding both fields prevents a future replacement from
    /// being mistaken for an already-flushed generation.
    flushed_idempotency: BTreeMap<(u64, String), (u64, SubmissionReceipt)>,
    /// Prefix length of `events` durably appended to the owning ledger.
    flushed_events: usize,
    /// Exact owning ledger sink. File ledgers bind by their opened path;
    /// independent `:memory:` handles additionally bind by handle identity.
    flushed_ledger: Option<String>,
    next_submission_ordinal: u64,
    next_ordinal: i64,
    next_flush_generation: i64,
}

/// The governor. `Send + Sync`: hot paths are mutex-guarded in-memory
/// state; ledger persistence is the explicit single-threaded flush.
pub struct Governor {
    inner: Mutex<Inner>,
    idle: Condvar,
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
            inner: Mutex::new(Inner::default()),
            idle: Condvar::new(),
        }
    }

    fn register_session(
        &self,
        token: CapabilityToken,
        gate: Option<Arc<CancelGate>>,
    ) -> Result<(), SessionError> {
        validate_resource("core-seconds grant", token.core_s)?;
        validate_resource("wall-seconds grant", token.wall_s)?;
        validate_resource("concurrent-cores grant", token.cores)?;
        let session = token.session.0;
        let mut g = self.inner.lock().expect("governor lock");
        if g.tokens.contains_key(&session) {
            return Err(SessionError::SessionAlreadyOpen { id: session });
        }
        g.meters.insert(session, SessionMeters::default());
        g.tokens.insert(session, token);
        if let Some(gate) = gate {
            g.gates.insert(session, gate);
        }
        Ok(())
    }

    /// Register a session's token (issuance). Session ids are single-use for
    /// the lifetime of this governor; duplicate registration fails closed.
    ///
    /// # Errors
    /// [`SessionError::InvalidResource`] when a floating-point grant is not
    /// finite and non-negative, or [`SessionError::SessionAlreadyOpen`] when
    /// the id is already registered. Rejection happens before any session
    /// state is mutated.
    pub fn open_session(&self, token: CapabilityToken) -> Result<(), SessionError> {
        self.register_session(token, None)
    }

    /// Register a session's token WITH its cancellation capability
    /// (bead gp3.13): the gate is owned by the governor from open, and
    /// level-3 memory pressure resolves it by `SessionId` — passing
    /// someone else's gate to a pressure action is unrepresentable.
    /// Sessions opened without a gate refuse level-3 pressure.
    ///
    /// # Errors
    /// [`SessionError::InvalidResource`] or
    /// [`SessionError::SessionAlreadyOpen`] as [`Governor::open_session`].
    pub fn open_session_gated(
        &self,
        token: CapabilityToken,
        gate: Arc<CancelGate>,
    ) -> Result<(), SessionError> {
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
    /// [`SessionError::UnknownSession`] or [`SessionError::InvalidResource`].
    pub fn charge(&self, session: SessionId, delta: Charge) -> Result<Enforcement, SessionError> {
        validate_resource("core-seconds charge", delta.core_s)?;
        validate_resource("wall-seconds charge", delta.wall_s)?;
        let mut g = self.inner.lock().expect("governor lock");
        let token = g
            .tokens
            .get(&session.0)
            .cloned()
            .ok_or(SessionError::UnknownSession { id: session.0 })?;
        let meters = g.meters.entry(session.0).or_default();
        let next_core_s = meters.core_s + delta.core_s;
        let next_wall_s = meters.wall_s + delta.wall_s;
        validate_resource("accumulated core-seconds", next_core_s)?;
        validate_resource("accumulated wall-seconds", next_wall_s)?;
        meters.core_s = next_core_s;
        meters.mem_peak_bytes = meters.mem_peak_bytes.max(delta.mem_peak_bytes);
        meters.wall_s = next_wall_s;
        #[allow(clippy::cast_precision_loss)]
        let checks: [(&'static str, f64, f64); 3] = [
            ("core-seconds", meters.core_s, token.core_s),
            (
                "memory-bytes",
                meters.mem_peak_bytes as f64,
                token.mem_bytes as f64,
            ),
            ("wall-seconds", meters.wall_s, token.wall_s),
        ];
        for (resource, used, granted) in checks {
            if used > granted * HARD_FACTOR {
                meters.paused += 1;
                return Ok(Enforcement::Paused {
                    resource,
                    used,
                    granted,
                    resume_hint: format!(
                        "checkpoint accepted; resume with a larger {resource} grant or a \
                         coarsened study — consumption and checkpoint are ledgered"
                    ),
                });
            }
        }
        for (resource, used, granted) in checks {
            if used >= granted {
                meters.throttled += 1;
                return Ok(Enforcement::Throttled {
                    resource,
                    used,
                    granted,
                });
            }
        }
        Ok(Enforcement::Ok)
    }

    /// Idempotency-keyed exactly-once execution: the first caller runs
    /// `work` and is charged; concurrent/repeat callers with the same key
    /// wait and receive `Duplicate` with the SAME receipt and NO charge.
    ///
    /// # Errors
    /// [`SessionError::UnknownSession`] for an unknown owner, or
    /// [`SessionError::Submission`] for a blank/oversized key or exhausted
    /// logical ordinal space.
    ///
    /// A panic in `work` is contained and committed as a terminal
    /// [`SubmitOutcome::Failed`] receipt. The same key never reruns implicitly:
    /// duplicates receive that same failure receipt and callers must choose a
    /// new idempotency key for an explicit retry.
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
        let scope = (session.0, idem_key.to_string());
        let ordinal = {
            let mut g = self.inner.lock().expect("governor lock");
            if !g.tokens.contains_key(&session.0) {
                return Err(SessionError::UnknownSession { id: session.0 });
            }
            loop {
                match g.idempotency.get(&scope) {
                    None => {
                        g.next_submission_ordinal = g
                            .next_submission_ordinal
                            .checked_add(1)
                            .ok_or_else(|| SessionError::Submission {
                                what: "submission ordinal space exhausted".to_string(),
                            })?;
                        let ordinal = g.next_submission_ordinal;
                        g.idempotency.insert(scope.clone(), IdemState::Pending);
                        break ordinal; // we own execution
                    }
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
                    Some(IdemState::Failed { receipt, what, .. }) => {
                        return Ok(SubmitOutcome::Failed {
                            receipt: *receipt,
                            what: what.clone(),
                        });
                    }
                    Some(IdemState::Pending) => {
                        g = self.idle.wait(g).expect("governor wait");
                    }
                }
            }
        };
        // Execute OUTSIDE the lock (work may be long). Catching here is
        // load-bearing: every Pending key must reach a terminal state and wake
        // its waiters even when caller-authored work unwinds.
        let completion = match std::panic::catch_unwind(std::panic::AssertUnwindSafe(work)) {
            Ok(charge) => match self.charge(session, charge) {
                Ok(enforcement) => SubmissionCompletion::Done(charge, enforcement),
                Err(error) => SubmissionCompletion::Failed(error.to_string()),
            },
            Err(payload) => SubmissionCompletion::Failed(panic_message(payload.as_ref())),
        };
        let receipt = submission_receipt(session, idem_key, &completion);
        let outcome;
        {
            let mut g = self.inner.lock().expect("governor lock");
            match completion {
                SubmissionCompletion::Done(charge, enforcement) => {
                    g.idempotency.insert(
                        scope,
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
                SubmissionCompletion::Failed(what) => {
                    g.idempotency.insert(
                        scope,
                        IdemState::Failed {
                            ordinal,
                            receipt,
                            what: what.clone(),
                        },
                    );
                    outcome = SubmitOutcome::Failed { receipt, what };
                }
            }
        }
        self.idle.notify_all();
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
    pub fn apply_memory_pressure(
        &self,
        session: SessionId,
        level: u8,
    ) -> Result<Vec<DegradationEvent>, SessionError> {
        if !(1..=3).contains(&level) {
            return Err(SessionError::InvalidPressureLevel { level });
        }
        let mut g = self.inner.lock().expect("governor lock");
        if !g.tokens.contains_key(&session.0) {
            return Err(SessionError::UnknownSession { id: session.0 });
        }
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
        let mut fired = Vec::new();
        for (i, step) in LADDER.iter().enumerate() {
            if i as u8 >= level {
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
                DegradationStep::PauseSerializeResume => {
                    gate.as_ref()
                        .expect("level-3 gate resolved above")
                        .request();
                    (
                        StepPhase::Requested,
                        "requested pause on the session-owned gate: solver checkpoints \
                         at the next tile boundary (SolverState snapshot to the ledger); \
                         complete only on acknowledge_pause with a checkpoint receipt"
                            .to_string(),
                    )
                }
            };
            g.next_ordinal += 1;
            let event = DegradationEvent {
                session,
                step: *step,
                pressure_level: level,
                phase,
                attribution,
                ordinal: g.next_ordinal,
            };
            if event.phase == StepPhase::Requested {
                g.pending_pause.insert(session.0, event.ordinal);
            }
            fired.push(event.clone());
            g.events.push(event);
        }
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
    pub fn acknowledge_pause(
        &self,
        session: SessionId,
        checkpoint_receipt: &str,
    ) -> Result<DegradationEvent, SessionError> {
        let mut g = self.inner.lock().expect("governor lock");
        if !g.tokens.contains_key(&session.0) {
            return Err(SessionError::UnknownSession { id: session.0 });
        }
        if checkpoint_receipt.trim().is_empty() {
            return Err(SessionError::Submission {
                what: "pause acknowledgement requires a non-empty checkpoint receipt".to_string(),
            });
        }
        let requested_ordinal = g
            .pending_pause
            .remove(&session.0)
            .ok_or(SessionError::NoPendingPause { id: session.0 })?;
        g.next_ordinal += 1;
        let event = DegradationEvent {
            session,
            step: DegradationStep::PauseSerializeResume,
            pressure_level: 3,
            phase: StepPhase::Complete,
            attribution: format!(
                "pause complete: checkpoint receipt {checkpoint_receipt:?} acknowledges \
                 the request at ordinal {requested_ordinal}"
            ),
            ordinal: g.next_ordinal,
        };
        g.events.push(event.clone());
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

    /// All recorded degradation events (deterministic ordinal order).
    #[must_use]
    pub fn events(&self) -> Vec<DegradationEvent> {
        self.inner.lock().expect("governor lock").events.clone()
    }

    /// Incrementally persist changed consumption snapshots plus new terminal
    /// idempotency and degradation events to the Design Ledger. A successful
    /// flush advances in-memory cursors, so repeating it without new governor
    /// state is a no-op. The batch is atomic and this governor has one owning
    /// ledger sink; cross-ledger replication belongs above this API.
    ///
    /// The method refuses to join an already-open ledger transaction: it could
    /// not know whether the caller later committed or rolled back, and marking
    /// the in-memory cursors in either case would lose data or duplicate it.
    /// (Single-threaded by design: fsqlite connections are `!Send`.)
    ///
    /// # Errors
    /// [`SessionError::Persistence`] wrapping the ledger error.
    pub fn flush_to_ledger(&self, ledger: &fs_ledger::Ledger) -> Result<(), SessionError> {
        let mut g = self.inner.lock().expect("governor lock");
        let sink_identity = ledger_sink_identity(ledger);
        if let Some(bound) = &g.flushed_ledger
            && bound != &sink_identity
        {
            return Err(SessionError::Persistence {
                what: format!(
                    "this governor is already bound to ledger sink {bound:?}; refusing split \
                     history into {sink_identity:?} and leaving every flush cursor unchanged"
                ),
            });
        }
        let flush_generation =
            g.next_flush_generation
                .checked_add(1)
                .ok_or_else(|| SessionError::Persistence {
                    what: "session flush generation space exhausted".to_string(),
                })?;
        let mut buffered = Vec::new();
        let mut meter_marks = Vec::new();
        for (id, m) in &g.meters {
            if g.flushed_meters
                .get(id)
                .is_some_and(|flushed| same_meter_snapshot(flushed, m))
            {
                continue;
            }
            let payload = format!(
                "{{\"schema\":\"fs-session-consumption-v2\",\"flush_generation\":{flush_generation},\"core_s\":{},\"mem_peak\":{},\"wall_s\":{},\"throttled\":{},\"paused\":{}}}",
                m.core_s, m.mem_peak_bytes, m.wall_s, m.throttled, m.paused,
            );
            buffered.push(BufferedLedgerEvent {
                session: id.to_be_bytes(),
                t: flush_generation,
                kind: "session.consumption",
                payload,
            });
            meter_marks.push((*id, m.clone()));
        }
        let mut idempotency_marks = Vec::new();
        for ((session, key), state) in &g.idempotency {
            let (ordinal, receipt, kind, payload) = match state {
                IdemState::Pending => continue,
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
                        "{{\"schema\":\"fs-session-idempotency-v2\",\"session\":{session},\"key\":\"{}\",\"receipt\":\"{receipt}\",\"core_s_bits\":\"{:016x}\",\"mem_peak_bytes\":{},\"wall_s_bits\":\"{:016x}\",\"enforcement\":{}}}",
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
                    what,
                } => (
                    *ordinal,
                    *receipt,
                    "session.idempotent-failure",
                    format!(
                        "{{\"schema\":\"fs-session-idempotency-v2\",\"session\":{session},\"key\":\"{}\",\"receipt\":\"{receipt}\",\"error\":\"{}\"}}",
                        json_escape(key),
                        json_escape(what),
                    ),
                ),
            };
            let generation = (ordinal, receipt);
            let scope = (*session, key.clone());
            if g.flushed_idempotency.get(&scope) == Some(&generation) {
                continue;
            }
            buffered.push(BufferedLedgerEvent {
                session: session.to_be_bytes(),
                t: i64::try_from(ordinal).map_err(|_| SessionError::Persistence {
                    what: format!("idempotency ordinal {ordinal} exceeds ledger i64"),
                })?,
                kind,
                payload,
            });
            idempotency_marks.push((scope, generation));
        }
        let flushed_events = g.flushed_events;
        let event_target = g.events.len();
        let new_events = g.events.get(flushed_events..).ok_or_else(|| {
            SessionError::Persistence {
                what: format!(
                    "session flush cursor {flushed_events} exceeds degradation event count {event_target}"
                ),
            }
        })?;
        for ev in new_events {
            let payload = format!(
                "{{\"step\":\"{:?}\",\"level\":{},\"phase\":\"{:?}\",\"attribution\":\"{}\"}}",
                ev.step,
                ev.pressure_level,
                ev.phase,
                json_escape(&ev.attribution),
            );
            buffered.push(BufferedLedgerEvent {
                session: ev.session.0.to_be_bytes(),
                t: ev.ordinal,
                kind: "session.degradation",
                payload,
            });
        }
        if buffered.is_empty() {
            return Ok(());
        }
        if ledger.in_transaction() {
            return Err(SessionError::Persistence {
                what: "session flush requires ownership of its atomic ledger transaction; an \
                       explicit transaction is already open and every flush cursor remains dirty"
                    .to_string(),
            });
        }
        let rows: Vec<_> = buffered.iter().map(BufferedLedgerEvent::as_row).collect();
        ledger
            .append_events(&rows)
            .map_err(|e| SessionError::Persistence {
                what: format!(
                    "atomic session event batch failed; every flush cursor remains dirty: {e}"
                ),
            })?;

        g.next_flush_generation = flush_generation;
        g.flushed_ledger = Some(sink_identity);
        for (id, meters) in meter_marks {
            g.flushed_meters.insert(id, meters);
        }
        for (scope, generation) in idempotency_marks {
            g.flushed_idempotency.insert(scope, generation);
        }
        g.flushed_events = event_target;
        Ok(())
    }
}

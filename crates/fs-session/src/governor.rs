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

/// Outcome of an idempotent submission.
#[derive(Debug, Clone, PartialEq)]
pub enum SubmitOutcome {
    /// This call executed the work.
    Executed {
        /// The charge recorded.
        charge: Charge,
        /// Execution ordinal (1 = first ever for this key).
        receipt: u64,
    },
    /// The key had already executed (or raced and lost): same receipt,
    /// NO additional charge.
    Duplicate {
        /// The original execution's receipt.
        receipt: u64,
    },
    /// The one attempted execution failed before a charge could be committed.
    /// The key remains terminal: all duplicates receive this same receipt and
    /// diagnosis, and an explicit retry requires a new key.
    Failed {
        /// The failed execution's receipt.
        receipt: u64,
        /// Panic payload or structured validation diagnosis.
        what: String,
    },
    /// Rejected with guidance before execution.
    Refused(Box<Guidance>),
}

#[derive(Debug, Default)]
struct SessionMeters {
    core_s: f64,
    mem_peak_bytes: u64,
    wall_s: f64,
    throttled: u32,
    paused: u32,
}

#[derive(Debug)]
enum IdemState {
    Pending,
    Done { receipt: u64, charge: Charge },
    Failed { receipt: u64, what: String },
}

enum SubmissionCompletion {
    Done(Charge),
    Failed(String),
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
    idempotency: BTreeMap<String, IdemState>,
    events: Vec<DegradationEvent>,
    next_receipt: u64,
    next_ordinal: i64,
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

    /// Register a session's token (issuance).
    ///
    /// # Errors
    /// [`SessionError::InvalidResource`] when a floating-point grant is not
    /// finite and non-negative. Rejection happens before any session state is
    /// mutated.
    pub fn open_session(&self, token: CapabilityToken) -> Result<(), SessionError> {
        validate_resource("core-seconds grant", token.core_s)?;
        validate_resource("wall-seconds grant", token.wall_s)?;
        validate_resource("concurrent-cores grant", token.cores)?;
        let mut g = self.inner.lock().expect("governor lock");
        g.meters.entry(token.session.0).or_default();
        g.tokens.insert(token.session.0, token);
        Ok(())
    }

    /// Register a session's token WITH its cancellation capability
    /// (bead gp3.13): the gate is owned by the governor from open, and
    /// level-3 memory pressure resolves it by `SessionId` — passing
    /// someone else's gate to a pressure action is unrepresentable.
    /// Sessions opened without a gate refuse level-3 pressure.
    ///
    /// # Errors
    /// [`SessionError::InvalidResource`] as [`Governor::open_session`].
    pub fn open_session_gated(
        &self,
        token: CapabilityToken,
        gate: Arc<CancelGate>,
    ) -> Result<(), SessionError> {
        let session = token.session.0;
        self.open_session(token)?;
        let mut g = self.inner.lock().expect("governor lock");
        g.gates.insert(session, gate);
        Ok(())
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
    /// [`SessionError::UnknownSession`].
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
    /// [`SessionError::UnknownSession`].
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
        {
            let mut g = self.inner.lock().expect("governor lock");
            if !g.tokens.contains_key(&session.0) {
                return Err(SessionError::UnknownSession { id: session.0 });
            }
            loop {
                match g.idempotency.get(idem_key) {
                    None => {
                        g.idempotency
                            .insert(idem_key.to_string(), IdemState::Pending);
                        break; // we own execution
                    }
                    Some(IdemState::Done { receipt, .. }) => {
                        return Ok(SubmitOutcome::Duplicate { receipt: *receipt });
                    }
                    Some(IdemState::Failed { receipt, what }) => {
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
        }
        // Execute OUTSIDE the lock (work may be long). Catching here is
        // load-bearing: every Pending key must reach a terminal state and wake
        // its waiters even when caller-authored work unwinds.
        let completion = match std::panic::catch_unwind(std::panic::AssertUnwindSafe(work)) {
            Ok(charge) => match self.charge(session, charge) {
                Ok(_) => SubmissionCompletion::Done(charge),
                Err(error) => SubmissionCompletion::Failed(error.to_string()),
            },
            Err(payload) => SubmissionCompletion::Failed(panic_message(payload.as_ref())),
        };
        let receipt;
        let outcome;
        {
            let mut g = self.inner.lock().expect("governor lock");
            g.next_receipt += 1;
            receipt = g.next_receipt;
            match completion {
                SubmissionCompletion::Done(charge) => {
                    g.idempotency
                        .insert(idem_key.to_string(), IdemState::Done { receipt, charge });
                    outcome = SubmitOutcome::Executed { charge, receipt };
                }
                SubmissionCompletion::Failed(what) => {
                    g.idempotency.insert(
                        idem_key.to_string(),
                        IdemState::Failed {
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

    /// The canonical idempotency key: agent-supplied key + content hash of
    /// the submitted program text.
    #[must_use]
    pub fn idempotency_key(agent_key: &str, program_text: &str) -> String {
        format!(
            "{agent_key}:{:016x}",
            fs_obs::fnv1a64(program_text.as_bytes())
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

    /// Persist consumption + degradation events to the Design Ledger
    /// (single-threaded by design: fsqlite connections are `!Send`).
    ///
    /// # Errors
    /// [`SessionError::Persistence`] wrapping the ledger error.
    pub fn flush_to_ledger(&self, ledger: &fs_ledger::Ledger) -> Result<(), SessionError> {
        let g = self.inner.lock().expect("governor lock");
        let persist = |e: &SessionError| e.to_string();
        for (id, m) in &g.meters {
            let payload = format!(
                "{{\"core_s\":{},\"mem_peak\":{},\"wall_s\":{},\"throttled\":{},\"paused\":{}}}",
                m.core_s, m.mem_peak_bytes, m.wall_s, m.throttled, m.paused
            );
            let session_bytes = id.to_be_bytes();
            ledger
                .append_event(&fs_ledger::EventRow {
                    session: Some(&session_bytes),
                    t: 0,
                    kind: "session.consumption",
                    payload: Some(&payload),
                })
                .map_err(|e| SessionError::Persistence {
                    what: format!("consumption event: {e}"),
                })?;
        }
        for (key, state) in &g.idempotency {
            let (receipt, kind, payload) = match state {
                IdemState::Pending => continue,
                IdemState::Done { receipt, charge } => (
                    *receipt,
                    "session.idempotent-execution",
                    format!(
                        "{{\"receipt\":{},\"core_s\":{},\"wall_s\":{}}}",
                        receipt, charge.core_s, charge.wall_s
                    ),
                ),
                IdemState::Failed { receipt, what } => (
                    *receipt,
                    "session.idempotent-failure",
                    format!("{{\"receipt\":{receipt},\"error\":{what:?}}}"),
                ),
            };
            ledger
                .append_event(&fs_ledger::EventRow {
                    session: None,
                    t: i64::try_from(receipt).unwrap_or(i64::MAX),
                    kind,
                    payload: Some(&payload),
                })
                .map_err(|e| SessionError::Persistence {
                    what: format!("idempotency event for {key}: {e}"),
                })?;
        }
        for ev in &g.events {
            let payload = format!(
                "{{\"step\":\"{:?}\",\"level\":{},\"phase\":\"{:?}\",\"attribution\":{:?}}}",
                ev.step, ev.pressure_level, ev.phase, ev.attribution
            );
            let session_bytes = ev.session.0.to_be_bytes();
            ledger
                .append_event(&fs_ledger::EventRow {
                    session: Some(&session_bytes),
                    t: ev.ordinal,
                    kind: "session.degradation",
                    payload: Some(&payload),
                })
                .map_err(|e| SessionError::Persistence {
                    what: persist(&SessionError::Persistence {
                        what: e.to_string(),
                    }),
                })?;
        }
        Ok(())
    }
}

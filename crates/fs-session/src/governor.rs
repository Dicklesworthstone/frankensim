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
use std::sync::{Condvar, Mutex};

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

/// A ledgered degradation event.
#[derive(Debug, Clone, PartialEq)]
pub struct DegradationEvent {
    /// The affected session.
    pub session: SessionId,
    /// Which ladder step fired.
    pub step: DegradationStep,
    /// Pressure level (1..=3) that triggered it.
    pub pressure_level: u8,
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
}

#[derive(Default)]
struct Inner {
    tokens: BTreeMap<u64, CapabilityToken>,
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
    pub fn open_session(&self, token: CapabilityToken) {
        let mut g = self.inner.lock().expect("governor lock");
        g.meters.entry(token.session.0).or_default();
        g.tokens.insert(token.session.0, token);
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
        let mut g = self.inner.lock().expect("governor lock");
        let token = g
            .tokens
            .get(&session.0)
            .cloned()
            .ok_or(SessionError::UnknownSession { id: session.0 })?;
        let meters = g.meters.entry(session.0).or_default();
        meters.core_s += delta.core_s;
        meters.mem_peak_bytes = meters.mem_peak_bytes.max(delta.mem_peak_bytes);
        meters.wall_s += delta.wall_s;
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
            if used > granted {
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
    /// # Panics
    /// If the governor mutex is poisoned (a prior panic in `work`).
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
                    Some(IdemState::Pending) => {
                        g = self.idle.wait(g).expect("governor wait");
                    }
                }
            }
        }
        // Execute OUTSIDE the lock (work may be long).
        let charge = work();
        let receipt;
        {
            let mut g = self.inner.lock().expect("governor lock");
            g.next_receipt += 1;
            receipt = g.next_receipt;
            g.idempotency
                .insert(idem_key.to_string(), IdemState::Done { receipt, charge });
        }
        self.idle.notify_all();
        // One charge, exactly once.
        let _ = self.charge(session, charge)?;
        Ok(SubmitOutcome::Executed { charge, receipt })
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

    /// Apply memory pressure at `level` (1..=3): ladder steps `1..=level`
    /// fire IN THE DECLARED ORDER, each recorded with attribution. The
    /// `PauseSerializeResume` step requests cancellation on the session's
    /// gate — the solver checkpoints at its next tile boundary (P7).
    ///
    /// # Errors
    /// [`SessionError::UnknownSession`].
    pub fn apply_memory_pressure(
        &self,
        session: SessionId,
        level: u8,
        gate: Option<&CancelGate>,
    ) -> Result<Vec<DegradationEvent>, SessionError> {
        let mut g = self.inner.lock().expect("governor lock");
        if !g.tokens.contains_key(&session.0) {
            return Err(SessionError::UnknownSession { id: session.0 });
        }
        let mut fired = Vec::new();
        for (i, step) in LADDER.iter().enumerate() {
            if i as u8 >= level {
                break;
            }
            let attribution = match step {
                DegradationStep::SpillColdArenas => {
                    "spilled coldest arenas (least-recently-touched first)".to_string()
                }
                DegradationStep::CoarsenAdaptively => {
                    "coarsened adaptive resolutions outside protected bands".to_string()
                }
                DegradationStep::PauseSerializeResume => {
                    if let Some(gate) = gate {
                        gate.request();
                    }
                    "requested pause: solver checkpoints at the next tile boundary \
                     (SolverState snapshot to the ledger)"
                        .to_string()
                }
            };
            g.next_ordinal += 1;
            let event = DegradationEvent {
                session,
                step: *step,
                pressure_level: level,
                attribution,
                ordinal: g.next_ordinal,
            };
            fired.push(event.clone());
            g.events.push(event);
        }
        Ok(fired)
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
            if let IdemState::Done { receipt, charge } = state {
                let payload = format!(
                    "{{\"receipt\":{},\"core_s\":{},\"wall_s\":{}}}",
                    receipt, charge.core_s, charge.wall_s
                );
                ledger
                    .append_event(&fs_ledger::EventRow {
                        session: None,
                        t: i64::try_from(*receipt).unwrap_or(i64::MAX),
                        kind: "session.idempotent-execution",
                        payload: Some(&payload),
                    })
                    .map_err(|e| SessionError::Persistence {
                        what: format!("idempotency event for {key}: {e}"),
                    })?;
            }
        }
        for ev in &g.events {
            let payload = format!(
                "{{\"step\":\"{:?}\",\"level\":{},\"attribution\":{:?}}}",
                ev.step, ev.pressure_level, ev.attribution
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

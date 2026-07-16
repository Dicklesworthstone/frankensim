//! fs-session conformance (the gp3.7 bead). Acceptance: budget
//! enforcement throttles/pauses with structured outcomes (never a silent
//! kill); double-submit with one idempotency key = one execution, one
//! charge (concurrency-stress-tested); estimate() accuracy tracked vs
//! actuals with a ledgered calibration report; the degradation ladder
//! fires in its declared order under synthetic memory pressure with
//! pause-serialize-resume equality; errors surface as ranked guidance.

use fs_exec::solver::{SolverState, codec};
use fs_exec::{CancelGate, DrainTracker, RunId};
use fs_plan::{CostModel, CostObservation, SealedCostModel};
use fs_session::{
    CalibrationReport, CapabilityToken, Charge, DegradationEvent, DegradationStep,
    DurableGovernorNonce, Enforcement, Estimate, Governor as SessionGovernor, Guidance,
    MAX_CAPABILITY_OP_BYTES, MAX_CAPABILITY_OPS, MAX_EVENT_PAGE_ROWS, MAX_FLUSH_ENCODED_BYTES,
    MAX_FLUSH_ROWS, MAX_IDEMPOTENCY_INPUT_BYTES, MAX_IDEMPOTENCY_KEYS_PER_SESSION,
    MAX_LEDGER_SCOPE_BYTES, MAX_SESSIONS_PER_GOVERNOR, MAX_SESSIONS_PER_SCOPE, ScopeFlushPermit,
    SessionError, SessionId, StepPhase, SubmitOutcome, estimate,
};
use std::collections::BTreeMap;
use std::sync::Arc;
use std::sync::atomic::{AtomicU32, AtomicU64, Ordering};

struct Governor {
    inner: SessionGovernor,
    next_mutation: AtomicU64,
}

impl Governor {
    fn new() -> Self {
        Self {
            inner: SessionGovernor::new(),
            next_mutation: AtomicU64::new(1),
        }
    }

    /// Durable facade twin (bead e61t6 fork (a)): flush-bound fixtures
    /// must run the durable submission protocol, since the ledger's
    /// preclaim doctrine refuses submission terminals without a
    /// pre-execution claim and positive permit.
    fn new_durable(ledger: &fs_ledger::Ledger, nonce: DurableGovernorNonce) -> Self {
        Self {
            inner: SessionGovernor::new_durable(ledger, nonce)
                .expect("durable conformance governor"),
            next_mutation: AtomicU64::new(1),
        }
    }

    fn next_key(&self, kind: &str, session: SessionId) -> String {
        let ordinal = self.next_mutation.fetch_add(1, Ordering::Relaxed);
        format!("legacy-conformance-{kind}-{}-{ordinal}", session.0)
    }

    fn idempotency_key(agent_key: &str, program_text: &str) -> Result<String, SessionError> {
        SessionGovernor::idempotency_key(agent_key, program_text)
    }

    fn open_session(&self, token: CapabilityToken) -> Result<ScopeFlushPermit, SessionError> {
        let open_id = self
            .inner
            .session_open_id(token.session, &self.next_key("open", token.session))?;
        self.inner
            .open_session(open_id, token)
            .map(|receipt| receipt.flush_permit())
    }

    fn open_session_gated(
        &self,
        token: CapabilityToken,
        gate: Arc<CancelGate>,
    ) -> Result<ScopeFlushPermit, SessionError> {
        let open_id = self
            .inner
            .session_open_id(token.session, &self.next_key("open", token.session))?;
        self.inner
            .open_session_gated(open_id, token, gate)
            .map(|receipt| receipt.flush_permit())
    }

    fn charge(&self, session: SessionId, delta: Charge) -> Result<Enforcement, SessionError> {
        let report_id = self
            .inner
            .meter_report_id(session, &self.next_key("meter", session))?;
        self.inner
            .charge(report_id, delta)
            .map(|receipt| receipt.enforcement().clone())
    }

    fn submit_once(
        &self,
        session: SessionId,
        key: &str,
        work: impl FnOnce() -> Charge,
    ) -> Result<SubmitOutcome, SessionError> {
        let request_id = self.inner.submission_request_id(session, key, key)?;
        self.inner.submit_once(request_id, work)
    }

    /// Durable-protocol submission (bead e61t6 fork (a)): files the
    /// pre-execution Pending claim + positive permit the ledger's
    /// preclaim doctrine (3a25a0d) demands, so the flushed terminal is
    /// admissible. The `key` doubles as the canonical program text —
    /// it must match the request-id digest input exactly.
    fn submit_once_durable(
        &self,
        ledger: &fs_ledger::Ledger,
        session: SessionId,
        key: &str,
        work: impl FnOnce() -> Charge,
    ) -> Result<SubmitOutcome, SessionError> {
        let request_id = self.inner.submission_request_id(session, key, key)?;
        self.inner
            .submit_once_durable(ledger, request_id, key, work)
    }

    fn apply_memory_pressure(
        &self,
        session: SessionId,
        level: u8,
    ) -> Result<Vec<DegradationEvent>, SessionError> {
        let action_id = self
            .inner
            .pressure_action_id(session, &self.next_key("pressure", session))?;
        self.inner
            .apply_memory_pressure(action_id, level)
            .map(|receipt| receipt.events().to_vec())
    }
}

impl core::ops::Deref for Governor {
    type Target = SessionGovernor;

    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-session/conformance\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn token(id: u64, core_s: f64, wall_s: f64) -> CapabilityToken {
    CapabilityToken {
        session: SessionId(id),
        ops: vec![
            "flux.*".to_string(),
            "ascent.*".to_string(),
            "xform.*".to_string(),
        ],
        core_s,
        mem_bytes: 64 * 1024 * 1024 * 1024,
        wall_s,
        cores: 16,
        ledger_scope: "main".to_string(),
    }
}

fn token_in_scope(id: u64, ledger_scope: &str) -> CapabilityToken {
    let mut token = token(id, 1.0e9, 1.0e9);
    token.ledger_scope = ledger_scope.to_string();
    token
}

fn durable_ledger_path(case: &str) -> String {
    static NEXT: AtomicU64 = AtomicU64::new(1);
    let ordinal = NEXT.fetch_add(1, Ordering::Relaxed);
    std::env::temp_dir()
        .join(format!(
            "fs-session-conformance-ujhp-{}-{ordinal}-{case}.ledger",
            std::process::id()
        ))
        .to_string_lossy()
        .into_owned()
}

fn solver_checkpoint(
    ledger: &fs_ledger::Ledger,
    request: fs_session::PauseRequestId,
    label: &str,
) -> fs_ledger::session_registry::SolverCheckpointReceipt {
    let authority = request.checkpoint_authority();
    let mut run_bytes = [0_u8; 8];
    run_bytes.copy_from_slice(&authority.as_bytes()[..8]);
    let run = RunId(u64::from_le_bytes(run_bytes));
    let gate = CancelGate::new_clock_free();
    let tracker = DrainTracker::new(run, &gate);
    let worker = tracker.register_worker().expect("fixture worker");
    gate.request();
    drop(worker);
    let report = tracker.finalize().expect("fixture run drained");
    let snapshot = fs_exec::solver::envelope::seal(0x4653_434b_5054, 1, run.0, label.as_bytes());
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

fn meter_snapshot_tuple(snapshot: fs_session::MeterSnapshot) -> (f64, u64, f64, u32, u32) {
    (
        snapshot.core_s,
        snapshot.mem_peak_bytes,
        snapshot.wall_s,
        snapshot.throttled,
        snapshot.paused,
    )
}

const SPOUT: &str = r#"(study "spout-laminar-v3"
  (seed 0x5EED0001) (versions (constellation :lock "2026-07"))
  (budget (wall 2h) (mem 96GiB) (qoi-rel-error 2e-2))
  (let lever (xform.level-set-velocity vessel :band 12mm :dof 4096))
  (ascent.optimize J :over lever :method (lbfgs :m 17)))"#;

fn lbm_cost_model() -> SealedCostModel {
    let obs: Vec<CostObservation> = (1..=12)
        .map(|k| CostObservation {
            size: f64::from(k) * 512.0,
            cost_s: 0.1 * f64::from(k) * 512.0,
        })
        .collect();
    SealedCostModel::provisional_unaudited(CostModel::fit(&obs).expect("fits"), "session-test")
}

#[test]
fn ss_001_token_bridges_into_static_admission() {
    let t = token(1, 3600.0, 7200.0);
    assert!(t.grants_op("flux.free-surface-lbm"));
    assert!(!t.grants_op("quantum.anneal"));
    assert!(!t.grants_op("flux."));
    assert!(!t.grants_op("flux..solve"));
    assert!(!t.grants_op("flux.*"));
    let cap = t.to_admission().expect("bounded grants project");
    assert!((cap.wall_s - 7200.0).abs() < f64::EPSILON);
    assert_eq!(cap.mem_bytes, t.mem_bytes);
    assert_eq!(cap.cores, t.cores);

    let mut boundary = token(2, 1.0, 1.0);
    boundary.mem_bytes = u64::MAX;
    boundary.cores = 9_007_199_254_740_993;
    let boundary_cap = boundary.to_admission().expect("bounded grants project");
    assert_eq!(boundary_cap.mem_bytes, u64::MAX);
    assert_eq!(boundary_cap.cores, 9_007_199_254_740_993);

    let mut mixed_invalid = t.clone();
    mixed_invalid.ops.push("*".to_string());
    assert!(!mixed_invalid.grants_op("flux.free-surface-lbm"));
    assert!(matches!(
        mixed_invalid.to_admission(),
        Err(SessionError::InvalidOperatorGrant { index: 3, .. })
    ));
    let mut over_count = t.clone();
    over_count.ops = (0..=MAX_CAPABILITY_OPS)
        .map(|index| format!("operator-{index}"))
        .collect();
    assert!(!over_count.grants_op("operator-0"));
    assert!(matches!(
        over_count.to_admission(),
        Err(SessionError::LimitExceeded {
            resource: "capability_operator_grants",
            ..
        })
    ));
    // The bridge feeds fs-ir admission directly.
    let node = fs_ir::sexpr::parse(SPOUT).expect("parses");
    let cx = fs_ir::admission::AdmissionContext {
        cost_freshness: None,
        router: None,
        chart_requirements: Vec::new(),
        cost_models: BTreeMap::new(),
        capability: Some(cap),
        regime: None,
        regime_policy: fs_ir::admission::RegimePolicy::Warn,
    };
    let report = fs_ir::admission::admit(&node, &cx);
    assert!(report.admitted, "{}", report.diagnosis());
    verdict(
        "ss-001",
        "token globs + admission bridge verified end to end",
    );
}

#[test]
fn ss_002_enforcement_throttles_then_pauses_never_kills() {
    let gov = Governor::new();
    gov.open_session(token(7, 100.0, 1e9)).expect("valid token");
    // Under the grant: Ok.
    let e1 = gov
        .charge(
            SessionId(7),
            Charge {
                core_s: 60.0,
                ..Charge::default()
            },
        )
        .expect("session");
    assert_eq!(e1, Enforcement::Ok);
    // Over the grant: Throttled (structured; work continues).
    let e2 = gov
        .charge(
            SessionId(7),
            Charge {
                core_s: 50.0,
                ..Charge::default()
            },
        )
        .expect("session");
    match e2 {
        Enforcement::Throttled {
            resource,
            used,
            granted,
        } => {
            assert_eq!(resource, "core-seconds");
            assert!(used > granted);
        }
        other => panic!("expected Throttled, got {other:?}"),
    }
    // Past the hard bound: Paused with a teaching resume hint.
    let e3 = gov
        .charge(
            SessionId(7),
            Charge {
                core_s: 50.0,
                ..Charge::default()
            },
        )
        .expect("session");
    match e3 {
        Enforcement::Paused {
            resource,
            resume_hint,
            ..
        } => {
            assert_eq!(resource, "core-seconds");
            assert!(
                resume_hint.contains("resume"),
                "hint must teach: {resume_hint}"
            );
            assert!(
                resume_hint.contains("checkpoint required")
                    && !resume_hint.contains("checkpoint accepted"),
                "charge() must not claim the separate checkpoint acknowledgement occurred: {resume_hint}"
            );
        }
        other => panic!("expected Paused, got {other:?}"),
    }
    let (core_s, _, _, throttled, paused) = gov.consumption(SessionId(7)).expect("meters");
    assert!((core_s - 160.0).abs() < 1e-9);
    assert_eq!((throttled, paused), (1, 1));
    // Unknown sessions are structured errors.
    assert!(gov.charge(SessionId(99), Charge::default()).is_err());
    verdict(
        "ss-002",
        "Ok -> Throttled -> Paused ladder with meters; no silent kills",
    );
}

#[test]
fn ss_002a_memory_enforcement_is_exact_above_f64_integer_precision() {
    const GRANT: u64 = (1_u64 << 53) + 1;

    let gov = Governor::new();
    let mut below = token(70, f64::MAX, f64::MAX);
    below.mem_bytes = GRANT;
    gov.open_session(below).expect("valid exact-byte token");
    assert_eq!(
        gov.charge(
            SessionId(70),
            Charge {
                mem_peak_bytes: GRANT - 1,
                ..Charge::default()
            }
        )
        .expect("exact below-grant charge"),
        Enforcement::Ok,
        "adjacent u64 byte values must not collapse through f64"
    );
    assert!(matches!(
        gov.charge(
            SessionId(70),
            Charge {
                mem_peak_bytes: GRANT,
                ..Charge::default()
            }
        )
        .expect("exact at-grant charge"),
        Enforcement::Throttled {
            resource: "memory-bytes",
            ..
        }
    ));

    let hard_boundary =
        u64::try_from(u128::from(GRANT) * 6 / 5).expect("fixture hard boundary fits u64");
    let at_hard = Governor::new();
    let mut at_hard_token = token(72, f64::MAX, f64::MAX);
    at_hard_token.mem_bytes = GRANT;
    at_hard
        .open_session(at_hard_token)
        .expect("valid exact-byte token");
    assert!(matches!(
        at_hard
            .charge(
                SessionId(72),
                Charge {
                    mem_peak_bytes: hard_boundary,
                    ..Charge::default()
                }
            )
            .expect("hard-boundary charge"),
        Enforcement::Throttled {
            resource: "memory-bytes",
            ..
        }
    ));

    let past_hard = Governor::new();
    let mut past_hard_token = token(73, f64::MAX, f64::MAX);
    past_hard_token.mem_bytes = GRANT;
    past_hard
        .open_session(past_hard_token)
        .expect("valid exact-byte token");
    assert!(matches!(
        past_hard
            .charge(
                SessionId(73),
                Charge {
                    mem_peak_bytes: hard_boundary + 1,
                    ..Charge::default()
                }
            )
            .expect("past-hard-boundary charge"),
        Enforcement::Paused {
            resource: "memory-bytes",
            ..
        }
    ));
}

#[test]
fn ss_002b_duplicate_session_open_preserves_original_authority_and_state() {
    let gov = Governor::new();
    let original_gate = Arc::new(CancelGate::new());
    let replacement_gate = Arc::new(CancelGate::new());
    let original = token(71, 100.0, 1_000.0);
    gov.open_session_gated(original.clone(), Arc::clone(&original_gate))
        .expect("original session");
    let charge = Charge {
        core_s: 7.0,
        mem_peak_bytes: 11,
        wall_s: 3.0,
    };
    let executed = gov
        .submit_once(SessionId(71), "immutable-session", || charge)
        .expect("original submission");
    let (receipt, enforcement) = match executed {
        SubmitOutcome::Executed {
            receipt,
            enforcement,
            ..
        } => (receipt, enforcement),
        other => panic!("original key must execute, got {other:?}"),
    };
    let meters = gov.consumption(SessionId(71)).expect("original meters");

    let mut replacement = token(71, 1.0, 2.0);
    replacement.ops = vec!["unrelated.*".to_string()];
    replacement.mem_bytes = 1;
    replacement.cores = 1;
    replacement.ledger_scope = "replacement-scope".to_string();
    assert_eq!(
        gov.open_session(replacement.clone()),
        Err(SessionError::SessionAlreadyOpen { id: 71 })
    );
    assert_eq!(
        gov.open_session_gated(replacement, Arc::clone(&replacement_gate)),
        Err(SessionError::SessionAlreadyOpen { id: 71 })
    );

    assert_eq!(
        gov.token(SessionId(71)).expect("original token retained"),
        original,
        "duplicate registration must not mutate authority"
    );
    assert_eq!(
        gov.consumption(SessionId(71))
            .expect("original meters retained"),
        meters,
        "duplicate registration must not reset or alter meters"
    );
    assert!(matches!(
        gov.submit_once(SessionId(71), "immutable-session", || {
            panic!("duplicate registration must not erase the original idempotency state")
        })
        .expect("original terminal key retained"),
        SubmitOutcome::Duplicate {
            receipt: duplicate_receipt,
            enforcement: duplicate_enforcement,
            ..
        } if duplicate_receipt == receipt && duplicate_enforcement == enforcement
    ));

    gov.apply_memory_pressure(SessionId(71), 3)
        .expect("the original gate remains bound");
    assert!(original_gate.is_requested());
    assert!(
        !replacement_gate.is_requested(),
        "a rejected replacement gate must never acquire session authority"
    );
}

#[test]
fn ss_003_idempotency_races_execute_exactly_once() {
    let gov = Arc::new(Governor::new());
    gov.open_session(token(3, 1e9, 1e9)).expect("valid token");
    let executions = Arc::new(AtomicU32::new(0));
    let key = Governor::idempotency_key("agent-a", SPOUT).expect("bounded canonical key");
    let mut handles = Vec::new();
    for _ in 0..16 {
        let gov = Arc::clone(&gov);
        let executions = Arc::clone(&executions);
        let key = key.clone();
        handles.push(std::thread::spawn(move || {
            gov.submit_once(SessionId(3), &key, || {
                executions.fetch_add(1, Ordering::SeqCst);
                std::thread::sleep(std::time::Duration::from_millis(5));
                Charge {
                    core_s: 10.0,
                    ..Charge::default()
                }
            })
            .expect("session known")
        }));
    }
    let outcomes: Vec<SubmitOutcome> = handles
        .into_iter()
        .map(|h| h.join().expect("join"))
        .collect();
    assert_eq!(
        executions.load(Ordering::SeqCst),
        1,
        "exactly one execution"
    );
    let executed: Vec<_> = outcomes
        .iter()
        .filter(|o| matches!(o, SubmitOutcome::Executed { .. }))
        .collect();
    assert_eq!(executed.len(), 1, "exactly one Executed outcome");
    let receipt = match executed[0] {
        SubmitOutcome::Executed { receipt, .. } => *receipt,
        _ => unreachable!(),
    };
    for o in &outcomes {
        match o {
            SubmitOutcome::Executed { .. } | SubmitOutcome::InFlight => {}
            SubmitOutcome::Duplicate { receipt: r, .. } => {
                assert_eq!(
                    *r, receipt,
                    "terminal duplicates share the original receipt"
                );
            }
            other => panic!("race produced an impossible outcome: {other:?}"),
        }
    }
    assert!(matches!(
        gov.submit_once(SessionId(3), &key, || panic!("terminal duplicate ran"))
            .expect("terminal state remains queryable"),
        SubmitOutcome::Duplicate { receipt: r, .. } if r == receipt
    ));
    // ONE charge only.
    let (core_s, ..) = gov.consumption(SessionId(3)).expect("meters");
    assert!(
        (core_s - 10.0).abs() < 1e-9,
        "double-submit must not double-spend"
    );
    // A different key executes independently.
    let other = gov
        .submit_once(SessionId(3), "agent-a:other", || Charge {
            core_s: 5.0,
            ..Charge::default()
        })
        .expect("ok");
    assert!(matches!(other, SubmitOutcome::Executed { .. }));
    verdict(
        "ss-003",
        "16-thread race: one execution, one charge, non-blocking InFlight or shared terminal receipt",
    );
}

#[test]
fn ss_003a_same_thread_reentrant_duplicate_returns_in_flight() {
    let gov = Governor::new();
    gov.open_session(token(30, 1e9, 1e9))
        .expect("reentrant fixture session");
    let executions = AtomicU32::new(0);
    let key = "same-thread-reentrant";

    let outer = gov
        .submit_once(SessionId(30), key, || {
            executions.fetch_add(1, Ordering::SeqCst);
            assert_eq!(
                gov.submit_once(SessionId(30), key, || {
                    panic!("a reentrant duplicate must never execute")
                })
                .expect("reentrant lookup"),
                SubmitOutcome::InFlight,
                "Pending is observable without waiting on the owning thread"
            );
            Charge {
                core_s: 2.0,
                ..Charge::default()
            }
        })
        .expect("outer execution returns");
    assert!(matches!(outer, SubmitOutcome::Executed { .. }));
    assert_eq!(executions.load(Ordering::SeqCst), 1);
    assert_eq!(
        gov.consumption(SessionId(30)).expect("meters").0.to_bits(),
        2.0_f64.to_bits()
    );
}

#[test]
fn ss_003b_idempotency_is_session_scoped_and_content_addressed() {
    let gov = Governor::new();
    gov.open_session(token(31, 1e9, 1e9)).expect("session 31");
    gov.open_session(token(32, 1e9, 1e9)).expect("session 32");
    let key =
        Governor::idempotency_key("agent:alpha", "program:beta").expect("bounded canonical key");
    assert!(key.starts_with("fs-session-idem-v3:"));
    assert_eq!(key.len(), "fs-session-idem-v3:".len() + 64);
    assert_eq!(
        key,
        Governor::idempotency_key("agent:alpha", "program:beta")
            .expect("deterministic canonical key")
    );
    assert_ne!(
        key,
        Governor::idempotency_key("agent", "alpha:program:beta").expect("bounded canonical key"),
        "length framing must distinguish delimiter-equivalent inputs"
    );

    let first = gov
        .submit_once(SessionId(31), &key, || Charge {
            core_s: 3.0,
            ..Charge::default()
        })
        .expect("first session executes");
    let second = gov
        .submit_once(SessionId(32), &key, || Charge {
            core_s: 5.0,
            ..Charge::default()
        })
        .expect("second session executes independently");
    let first_receipt = match first {
        SubmitOutcome::Executed { receipt, .. } => receipt,
        other => panic!("first session must execute, got {other:?}"),
    };
    let second_receipt = match second {
        SubmitOutcome::Executed { receipt, .. } => receipt,
        other => panic!("second session must execute, got {other:?}"),
    };
    assert_ne!(
        first_receipt, second_receipt,
        "the owning session is part of receipt identity"
    );
    assert_eq!(
        gov.consumption(SessionId(31)).unwrap().0.to_bits(),
        3.0_f64.to_bits()
    );
    assert_eq!(
        gov.consumption(SessionId(32)).unwrap().0.to_bits(),
        5.0_f64.to_bits()
    );
    assert!(matches!(
        gov.submit_once(SessionId(31), &key, || panic!("duplicate ran"))
            .expect("duplicate returns"),
        SubmitOutcome::Duplicate { receipt, .. } if receipt == first_receipt
    ));
    assert!(
        gov.submit_once(SessionId(31), "   ", Charge::default)
            .is_err(),
        "blank idempotency keys must fail before execution"
    );
    let oversized = "x".repeat(MAX_IDEMPOTENCY_INPUT_BYTES + 1);
    assert!(matches!(
        Governor::idempotency_key("agent", &oversized),
        Err(SessionError::LimitExceeded {
            resource: "idempotency_program_text_bytes",
            ..
        })
    ));
}

#[test]
#[allow(clippy::too_many_lines)] // One receipt-identity conformance scenario.
fn ss_003c_receipts_bind_charge_failure_and_enforcement() {
    fn opened(core_s: f64) -> Governor {
        let governor = Governor::new();
        governor
            .open_session(token(33, core_s, 1e9))
            .expect("valid receipt-test token");
        governor
    }

    let key = "receipt-binding";
    let positive_zero_charge = Charge {
        core_s: 0.0,
        mem_peak_bytes: 7,
        wall_s: 1.0,
    };
    let positive_governor = opened(10.0);
    let positive_request = positive_governor
        .inner
        .submission_request_id(SessionId(33), key, key)
        .expect("typed request");
    let positive = positive_governor
        .submit_once(SessionId(33), key, || positive_zero_charge)
        .expect("positive-zero submit");
    let (positive_receipt, positive_admission, positive_meter, positive_enforcement) =
        match positive {
            SubmitOutcome::Executed {
                admission_ordinal,
                receipt,
                meter_receipt,
                enforcement,
                ..
            } => (receipt, admission_ordinal, meter_receipt, enforcement),
            other => panic!("expected execution, got {other:?}"),
        };
    assert_eq!(positive_enforcement, Enforcement::Ok);
    assert!(positive_receipt.matches_success(
        positive_request,
        "main",
        positive_admission,
        positive_zero_charge,
        &positive_meter,
    ));
    assert!(
        !positive_receipt.matches_success(
            positive_request,
            "other-scope",
            positive_admission,
            positive_zero_charge,
            &positive_meter,
        ),
        "receipt v3 binds the immutable ledger scope"
    );
    assert!(
        !positive_receipt.matches_success(
            positive_request,
            "main",
            positive_admission,
            Charge {
                mem_peak_bytes: positive_zero_charge.mem_peak_bytes + 1,
                ..positive_zero_charge
            },
            &positive_meter,
        ),
        "the exact u64 memory charge participates in receipt identity"
    );

    let negative_zero_charge = Charge {
        core_s: -0.0,
        ..positive_zero_charge
    };
    let negative = opened(10.0)
        .submit_once(SessionId(33), key, || negative_zero_charge)
        .expect("negative-zero submit");
    let negative_receipt = match negative {
        SubmitOutcome::Executed { receipt, .. } => receipt,
        other => panic!("expected execution, got {other:?}"),
    };
    assert_ne!(
        positive_receipt, negative_receipt,
        "bit-exact charge fields are receipt semantics"
    );

    let failed_governor = opened(10.0);
    let failed_request = failed_governor
        .inner
        .submission_request_id(SessionId(33), key, key)
        .expect("typed failure request");
    let failed = failed_governor
        .submit_once(SessionId(33), key, || panic!("receipt-bound failure"))
        .expect("panic is a terminal outcome");
    let failed_receipt = match failed {
        SubmitOutcome::Failed {
            admission_ordinal,
            receipt,
            evidence,
        } => {
            assert!(receipt.matches_failure(failed_request, "main", admission_ordinal, &evidence));
            assert!(evidence.preview().contains("receipt-bound failure"));
            receipt
        }
        other => panic!("expected failed receipt, got {other:?}"),
    };
    assert_ne!(positive_receipt, failed_receipt);

    // The flush-bound instance runs the DURABLE protocol (e61t6 fork
    // (a)): the ledger refuses submission terminals without a durable
    // pre-execution claim, so in-memory receipts stay in memory and
    // only claim-backed submissions flush.
    let throttled_path = durable_ledger_path("receipt-binding-throttled");
    let ledger = fs_ledger::Ledger::open(&throttled_path).expect("receipt ledger");
    let throttled = Governor::new_durable(&ledger, DurableGovernorNonce::from_bytes([0xC3; 32]));
    let throttled_permit = throttled
        .open_session(token(33, 1.0, 1e9))
        .expect("valid receipt-test token");
    // Durable prerequisite: the session-open terminal must be recorded
    // before durable submissions may run (RecoveryRequired otherwise).
    throttled
        .flush_scope_to_ledger(&throttled_permit, &ledger)
        .expect("open prerequisite terminal");
    let charge = Charge {
        core_s: 1.0,
        ..Charge::default()
    };
    let executed = throttled
        .submit_once_durable(&ledger, SessionId(33), key, || charge)
        .expect("throttled work still completes");
    let (receipt, admission_ordinal, meter_receipt, enforcement) = match executed {
        SubmitOutcome::Executed {
            admission_ordinal,
            receipt,
            meter_receipt,
            enforcement,
            ..
        } => (receipt, admission_ordinal, meter_receipt, enforcement),
        other => panic!("expected execution, got {other:?}"),
    };
    assert!(matches!(enforcement, Enforcement::Throttled { .. }));
    let throttled_request = throttled
        .inner
        .submission_request_id(SessionId(33), key, key)
        .expect("typed throttled request");
    assert!(receipt.matches_success(
        throttled_request,
        "main",
        admission_ordinal,
        charge,
        &meter_receipt
    ));
    assert!(matches!(
        throttled
            .submit_once_durable(&ledger, SessionId(33), key, || panic!("duplicate ran"))
            .expect("duplicate"),
        SubmitOutcome::Duplicate {
            receipt: duplicate_receipt,
            enforcement: Enforcement::Throttled { .. },
            ..
        } if duplicate_receipt == receipt
    ));

    throttled
        .flush_scope_to_ledger(&throttled_permit, &ledger)
        .expect("strict JSON receipt event");
    // Durable protocol: the flushed batch carries the pre-execution
    // claim's rows alongside the terminal receipt event, so the event
    // count is read from the observed durable layout rather than the
    // retired in-memory two-event shape.
    assert!(ledger.table_count("events").unwrap() >= 2);
    assert!(ledger.lint().unwrap().is_clean());
}

#[test]
fn ss_003d_ledger_flush_is_atomic_incremental_and_retryable() {
    // Durable protocol throughout (e61t6 fork (a)): the flush target is
    // the governor's recovery ledger and submissions carry pre-execution
    // claims, so flushed terminals satisfy the preclaim doctrine.
    let path = durable_ledger_path("flush-atomic");
    let ledger = fs_ledger::Ledger::open(&path).expect("flush ledger");
    let gov = Governor::new_durable(&ledger, DurableGovernorNonce::from_bytes([0xD4; 32]));
    let permit = gov
        .open_session(token(34, 1e9, 1e9))
        .expect("flush-test session");
    gov.flush_scope_to_ledger(&permit, &ledger)
        .expect("open prerequisite terminal");
    let after_open = ledger.table_count("events").unwrap();
    let key = "flush-once";
    assert!(matches!(
        gov.submit_once_durable(&ledger, SessionId(34), key, || Charge {
            core_s: 2.0,
            mem_peak_bytes: 5,
            wall_s: 1.0,
        })
        .expect("submission"),
        SubmitOutcome::Executed { .. }
    ));

    ledger.begin().expect("caller transaction");
    let refused = gov
        .flush_scope_to_ledger(&permit, &ledger)
        .expect_err("flush cannot promise durability inside a caller transaction");
    assert!(matches!(refused, SessionError::Persistence { .. }));
    assert_eq!(ledger.table_count("events").unwrap(), after_open);
    ledger.rollback().expect("caller rollback");

    gov.flush_scope_to_ledger(&permit, &ledger)
        .expect("the refused batch remains fully dirty for retry");
    let after_submission = ledger.table_count("events").unwrap();
    assert_eq!(
        after_submission,
        after_open + 1,
        "the claim-backed terminal submission receipt appends exactly once"
    );
    gov.flush_scope_to_ledger(&permit, &ledger)
        .expect("unchanged no-op flush");
    assert_eq!(
        ledger.table_count("events").unwrap(),
        after_submission,
        "repeated flush must not duplicate semantic events"
    );
    assert!(matches!(
        gov.submit_once_durable(&ledger, SessionId(34), key, || panic!("duplicate ran"))
            .expect("terminal duplicate"),
        SubmitOutcome::Duplicate { .. }
    ));
    gov.flush_scope_to_ledger(&permit, &ledger)
        .expect("duplicate observation is not a new event");
    assert_eq!(ledger.table_count("events").unwrap(), after_submission);

    gov.charge(
        SessionId(34),
        Charge {
            core_s: 3.0,
            ..Charge::default()
        },
    )
    .expect("new consumption");
    let foreign_ledger = fs_ledger::Ledger::open(":memory:").expect("foreign ledger");
    assert!(matches!(
        gov.flush_scope_to_ledger(&permit, &foreign_ledger),
        Err(SessionError::LedgerScopeSinkMismatch { .. })
    ));
    assert_eq!(
        foreign_ledger.table_count("events").unwrap(),
        0,
        "a governor must not split its event history across ledger sinks"
    );
    gov.flush_scope_to_ledger(&permit, &ledger)
        .expect("sink refusal leaves the changed meter dirty for its owning ledger");
    assert_eq!(ledger.table_count("events").unwrap(), after_submission + 1);
    gov.apply_memory_pressure(SessionId(34), 1)
        .expect("one degradation event");
    gov.flush_scope_to_ledger(&permit, &ledger)
        .expect("new degradation event appends once");
    assert_eq!(ledger.table_count("events").unwrap(), after_submission + 2);
    gov.flush_scope_to_ledger(&permit, &ledger)
        .expect("second no-op flush");
    assert_eq!(ledger.table_count("events").unwrap(), after_submission + 2);
    assert!(ledger.lint().unwrap().is_clean());
}

#[test]
#[allow(clippy::too_many_lines)] // Complete typed open/meter replay and conflict matrix.
fn ss_003m_typed_open_and_meter_authorities_replay_without_double_spend() {
    let governor = Arc::new(SessionGovernor::new());
    let session = SessionId(35);
    let capability = token(35, 100.0, 1e9);
    let open_id = governor
        .session_open_id(session, "typed-open-35")
        .expect("bounded open authority");
    let open_barrier = Arc::new(std::sync::Barrier::new(3));
    let mut open_workers = Vec::new();
    for _ in 0..2 {
        let governor = Arc::clone(&governor);
        let capability = capability.clone();
        let open_barrier = Arc::clone(&open_barrier);
        open_workers.push(std::thread::spawn(move || {
            open_barrier.wait();
            governor
                .open_session(open_id, capability)
                .expect("commit or replay open")
        }));
    }
    open_barrier.wait();
    let opened = open_workers.remove(0).join().expect("first open worker");
    let replayed = open_workers.remove(0).join().expect("second open worker");
    assert_eq!(opened, replayed);
    assert_eq!(
        governor
            .open_session(open_id, capability.clone())
            .expect("lost open response replays"),
        opened
    );
    let foreign = SessionGovernor::new();
    assert!(matches!(
        foreign.open_session(open_id, capability.clone()),
        Err(SessionError::MutationAuthorityMismatch {
            kind: "session-open",
            ..
        })
    ));
    let mut altered = capability;
    altered.core_s = 99.0;
    assert!(matches!(
        governor.open_session(open_id, altered),
        Err(SessionError::MutationConflict {
            kind: "session-open",
            ..
        })
    ));

    let report_id = governor
        .meter_report_id(session, "typed-meter-35")
        .expect("bounded report authority");
    let charge = Charge {
        core_s: 7.0,
        mem_peak_bytes: 11,
        wall_s: 3.0,
    };
    let barrier = Arc::new(std::sync::Barrier::new(3));
    let mut workers = Vec::new();
    for _ in 0..2 {
        let governor = Arc::clone(&governor);
        let barrier = Arc::clone(&barrier);
        workers.push(std::thread::spawn(move || {
            barrier.wait();
            governor
                .charge(report_id, charge)
                .expect("commit or exact replay")
        }));
    }
    barrier.wait();
    let first = workers.remove(0).join().expect("first meter worker");
    let second = workers.remove(0).join().expect("second meter worker");
    assert_eq!(first, second);
    assert_eq!(
        governor.consumption(session).expect("meter state"),
        (7.0, 11, 3.0, 0, 0),
        "concurrent duplicate changes the meter exactly once"
    );
    assert!(matches!(
        governor.charge(
            report_id,
            Charge {
                core_s: -0.0,
                ..charge
            }
        ),
        Err(SessionError::MutationConflict {
            kind: "meter-report",
            ..
        })
    ));
    assert!(matches!(
        foreign.charge(report_id, charge),
        Err(SessionError::MutationAuthorityMismatch {
            kind: "meter-report",
            ..
        })
    ));
    let zero_report = governor
        .meter_report_id(session, "signed-zero-report")
        .expect("signed-zero report authority");
    governor
        .charge(
            zero_report,
            Charge {
                core_s: 0.0,
                ..Charge::default()
            },
        )
        .expect("positive zero report");
    assert!(matches!(
        governor.charge(
            zero_report,
            Charge {
                core_s: -0.0,
                ..Charge::default()
            }
        ),
        Err(SessionError::MutationConflict {
            kind: "meter-report",
            ..
        })
    ));

    let gated = Arc::new(SessionGovernor::new());
    let gated_token = token(36, 100.0, 1e9);
    let gated_id = gated
        .session_open_id(SessionId(36), "typed-gated-open")
        .expect("gated open authority");
    let gate = Arc::new(CancelGate::new());
    let gated_barrier = Arc::new(std::sync::Barrier::new(3));
    let mut gated_workers = Vec::new();
    for _ in 0..2 {
        let gated = Arc::clone(&gated);
        let gated_token = gated_token.clone();
        let gate = Arc::clone(&gate);
        let gated_barrier = Arc::clone(&gated_barrier);
        gated_workers.push(std::thread::spawn(move || {
            gated_barrier.wait();
            gated
                .open_session_gated(gated_id, gated_token, gate)
                .expect("commit or replay gated open")
        }));
    }
    gated_barrier.wait();
    let gated_receipt = gated_workers.remove(0).join().expect("first gated worker");
    assert_eq!(
        gated_workers.remove(0).join().expect("second gated worker"),
        gated_receipt
    );
    let foreign_gated = SessionGovernor::new();
    assert!(matches!(
        foreign_gated.open_session_gated(gated_id, gated_token.clone(), Arc::clone(&gate),),
        Err(SessionError::MutationAuthorityMismatch {
            kind: "session-open",
            ..
        })
    ));
    assert!(matches!(
        gated.open_session_gated(gated_id, gated_token.clone(), Arc::new(CancelGate::new())),
        Err(SessionError::MutationConflict {
            kind: "session-open",
            ..
        })
    ));
    let draining_submission = gated
        .submission_request_id(SessionId(36), "draining-caller", "draining-program")
        .expect("submission authority before external cancellation");
    gate.request();
    assert_eq!(
        gated
            .open_session_gated(gated_id, gated_token, Arc::clone(&gate))
            .expect("exact replay precedes stale gate validation"),
        gated_receipt
    );
    assert!(matches!(
        gated.submit_once(draining_submission, || panic!("cancelled-gate work ran")),
        Err(SessionError::SessionGateDraining {
            id: 36,
            generation: 0,
        })
    ));
}

#[test]
fn ss_003n_pressure_action_replays_across_pause_lifecycle_and_refuses_stale_ids() {
    let governor = Arc::new(SessionGovernor::new());
    let session = SessionId(37);
    let capability = token(37, 1e9, 1e9);
    let open_id = governor
        .session_open_id(session, "pressure-open")
        .expect("open authority");
    let gate = Arc::new(CancelGate::new());
    let permit = governor
        .open_session_gated(open_id, capability, Arc::clone(&gate))
        .expect("gated open")
        .flush_permit();
    let action_id = governor
        .pressure_action_id(session, "pressure-action")
        .expect("action authority");
    let foreign = SessionGovernor::new();
    assert!(matches!(
        foreign.apply_memory_pressure(action_id, 3),
        Err(SessionError::MutationAuthorityMismatch {
            kind: "pressure-action",
            ..
        })
    ));
    let stale_unused = governor
        .pressure_action_id(session, "stale-unused")
        .expect("unused generation-zero authority");
    let stale_meter = governor
        .meter_report_id(session, "stale-unused-meter")
        .expect("unused generation-zero meter authority");
    let stale_submission = governor
        .submission_request_id(session, "stale-submission", "stale-program")
        .expect("unused generation-zero submission authority");
    let committed_meter = governor
        .meter_report_id(session, "committed-generation-zero-meter")
        .expect("committed generation-zero meter authority");
    let committed_meter_receipt = governor
        .charge(committed_meter, Charge::default())
        .expect("commit generation-zero meter");
    let barrier = Arc::new(std::sync::Barrier::new(3));
    let mut workers = Vec::new();
    for _ in 0..2 {
        let governor = Arc::clone(&governor);
        let barrier = Arc::clone(&barrier);
        workers.push(std::thread::spawn(move || {
            barrier.wait();
            governor
                .apply_memory_pressure(action_id, 3)
                .expect("commit or replay pressure")
        }));
    }
    barrier.wait();
    let first = workers.remove(0).join().expect("first pressure worker");
    let concurrent_replay = workers.remove(0).join().expect("second pressure worker");
    assert_eq!(first, concurrent_replay);
    assert_eq!(
        governor
            .apply_memory_pressure(action_id, 3)
            .expect("pending replay"),
        first
    );
    assert!(matches!(
        governor.apply_memory_pressure(action_id, 2),
        Err(SessionError::MutationConflict {
            kind: "pressure-action",
            ..
        })
    ));
    let request_id = first
        .events()
        .last()
        .and_then(|event| event.pause_request_id)
        .expect("pause request");
    let ledger = fs_ledger::Ledger::open(":memory:").expect("fixture ledger");
    let checkpoint = solver_checkpoint(&ledger, request_id, "typed-pressure-checkpoint");
    let acknowledgement = governor
        .acknowledge_pause(request_id, &ledger, &checkpoint)
        .expect("pause completion");
    assert_eq!(
        governor
            .apply_memory_pressure(action_id, 3)
            .expect("ready-to-resume replay"),
        first
    );
    governor
        .activate_resume(&acknowledgement)
        .expect("activate generation one");
    assert_eq!(
        governor
            .apply_memory_pressure(action_id, 3)
            .expect("activated replay"),
        first
    );
    assert!(matches!(
        governor.apply_memory_pressure(stale_unused, 1),
        Err(SessionError::StaleMutationGeneration {
            kind: "pressure-action",
            supplied: 0,
            current: 1,
            ..
        })
    ));
    assert!(matches!(
        governor.charge(stale_meter, Charge::default()),
        Err(SessionError::StaleMutationGeneration {
            kind: "meter-report",
            supplied: 0,
            current: 1,
            ..
        })
    ));
    assert_eq!(
        governor
            .charge(committed_meter, Charge::default())
            .expect("known old-generation report replays"),
        committed_meter_receipt
    );
    let executions = AtomicU32::new(0);
    assert!(matches!(
        governor.submit_once(stale_submission, || {
            executions.fetch_add(1, Ordering::SeqCst);
            Charge::default()
        }),
        Err(SessionError::StaleMutationGeneration {
            kind: "submission-request",
            supplied: 0,
            current: 1,
            ..
        })
    ));
    assert_eq!(executions.load(Ordering::SeqCst), 0);
    assert_eq!(
        governor
            .events_page(&permit, 0, MAX_EVENT_PAGE_ROWS)
            .expect("bounded events")
            .len(),
        4,
        "three action events plus one completion are retained exactly once"
    );

    let low_pressure = SessionGovernor::new();
    let low_session = SessionId(137);
    let low_open = low_pressure
        .session_open_id(low_session, "low-pressure-open")
        .expect("low-pressure open authority");
    low_pressure
        .open_session(low_open, token(137, 1e9, 1e9))
        .expect("ungated low-pressure session");
    for level in 1..=2 {
        let action = low_pressure
            .pressure_action_id(low_session, &format!("level-{level}"))
            .expect("low-pressure authority");
        let receipt = low_pressure
            .apply_memory_pressure(action, level)
            .expect("first low-pressure action");
        assert_eq!(
            low_pressure
                .apply_memory_pressure(action, level)
                .expect("exact low-pressure replay"),
            receipt
        );
    }
}

#[test]
#[allow(clippy::too_many_lines)] // One barrier-controlled causal inversion plus its durable replay proof.
fn ss_003o_submission_meter_commit_order_is_completion_order_and_atomic() {
    let path = durable_ledger_path("causal-inversion");
    let nonce = DurableGovernorNonce::from_bytes([0x38; 32]);
    let ledger = fs_ledger::Ledger::open(&path).expect("causal file ledger");
    let governor =
        Arc::new(SessionGovernor::new_durable(&ledger, nonce).expect("durable causal governor"));
    let session = SessionId(38);
    let session_token = token(38, 20.0, 1e9);
    let open_id = governor
        .session_open_id(session, "causal-open")
        .expect("open authority");
    let open_receipt = governor
        .open_session(open_id, session_token.clone())
        .expect("open session");
    let permit = open_receipt.flush_permit();
    let open_flush = governor
        .flush_scope_to_ledger(&permit, &ledger)
        .expect("durable open prerequisite");
    assert_eq!(open_flush.committed_terminals, 1);
    assert_eq!(open_flush.appended_rows, 1);
    assert!(!open_flush.remaining_dirty);
    let request_a = governor
        .submission_request_id(session, "caller-a", "program-a")
        .expect("request A");
    let conflicting_a = governor
        .submission_request_id(session, "caller-a", "different-program")
        .expect("conflicting request A");
    let request_b = governor
        .submission_request_id(session, "caller-b", "program-b")
        .expect("request B");
    let foreign = SessionGovernor::new();
    assert!(matches!(
        foreign.submit_once(request_a, Charge::default),
        Err(SessionError::MutationAuthorityMismatch {
            kind: "submission-request",
            ..
        })
    ));
    let (started_tx, started_rx) = std::sync::mpsc::channel();
    let (release_tx, release_rx) = std::sync::mpsc::channel();
    let worker = {
        let governor = Arc::clone(&governor);
        let worker_path = path.clone();
        std::thread::spawn(move || {
            let worker_ledger =
                fs_ledger::Ledger::open(&worker_path).expect("worker causal ledger handle");
            governor.submit_once_durable(&worker_ledger, request_a, "program-a", || {
                started_tx.send(()).expect("test receiver");
                release_rx.recv().expect("release sender");
                Charge {
                    core_s: 15.0,
                    ..Charge::default()
                }
            })
        })
    };
    started_rx.recv().expect("request A admitted first");
    assert!(matches!(
        governor.submit_once_durable(&ledger, conflicting_a, "different-program", Charge::default,),
        Err(SessionError::MutationConflict {
            kind: "submission",
            ..
        })
    ));
    let outcome_b = governor
        .submit_once_durable(&ledger, request_b, "program-b", || Charge {
            core_s: 5.0,
            ..Charge::default()
        })
        .expect("request B completes first");
    release_tx.send(()).expect("release request A");
    let outcome_a = worker
        .join()
        .expect("request A worker")
        .expect("request A completes second");
    let (admission_a, receipt_a, meter_a) = match outcome_a {
        SubmitOutcome::Executed {
            admission_ordinal,
            receipt,
            meter_receipt,
            ..
        } => (admission_ordinal, receipt, meter_receipt),
        other => panic!("request A must execute, got {other:?}"),
    };
    let (admission_b, receipt_b, meter_b) = match outcome_b {
        SubmitOutcome::Executed {
            admission_ordinal,
            receipt,
            meter_receipt,
            ..
        } => (admission_ordinal, receipt, meter_receipt),
        other => panic!("request B must execute, got {other:?}"),
    };
    assert_eq!((admission_a, admission_b), (1, 2));
    assert_eq!((meter_b.commit_ordinal(), meter_a.commit_ordinal()), (1, 2));
    assert!(admission_a < admission_b, "A was admitted first");
    assert!(
        meter_b.commit_ordinal() < meter_a.commit_ordinal(),
        "B completed and committed first"
    );
    assert_eq!(
        meter_b.after(),
        meter_a.before(),
        "causal receipts form an exact pre/post chain"
    );
    assert_eq!(meter_b.enforcement(), &Enforcement::Ok);
    assert!(matches!(
        meter_a.enforcement(),
        Enforcement::Throttled {
            resource: "core-seconds",
            used,
            granted,
        } if *used == 20.0 && *granted == 20.0
    ));
    let before_duplicate = governor.consumption(session).expect("meter state");
    assert!(matches!(
        governor
            .submit_once_durable(&ledger, request_a, "program-a", || {
                panic!("duplicate reran")
            })
            .expect("exact replay"),
        SubmitOutcome::Duplicate {
            enforcement,
            meter_receipt,
            ..
        } if meter_receipt == meter_a && enforcement == *meter_a.enforcement()
    ));
    assert_eq!(
        governor.consumption(session).expect("meter state"),
        before_duplicate
    );

    let terminal_flush = governor
        .flush_scope_to_ledger(&permit, &ledger)
        .expect("causal terminal flush");
    assert_eq!(terminal_flush.committed_terminals, 2);
    assert_eq!(terminal_flush.appended_rows, 2);
    assert!(!terminal_flush.remaining_dirty);
    assert_eq!(
        ledger.table_count("events").expect("event count"),
        3,
        "one open receipt and two self-contained causal terminal receipts"
    );
    assert!(ledger.lint().expect("ledger lint").is_clean());
    let durable_counts = (
        ledger.table_count("session_claims").unwrap(),
        ledger.table_count("session_terminals").unwrap(),
        ledger.table_count("session_terminal_events").unwrap(),
        ledger.table_count("session_flush_batches").unwrap(),
        ledger.table_count("session_flush_batch_members").unwrap(),
        ledger.table_count("events").unwrap(),
    );
    assert_eq!(durable_counts, (3, 3, 3, 2, 3, 3));
    drop(governor);
    drop(ledger);

    let ledger = fs_ledger::Ledger::open(&path).expect("reopened causal ledger");
    let governor = SessionGovernor::new_durable(&ledger, nonce).expect("reopened governor");
    assert!(matches!(
        governor.recover_open(&ledger, open_id, token(39, 20.0, 1e9), None),
        Err(SessionError::MutationAuthorityMismatch {
            kind: "session-open",
            ..
        })
    ));
    let mut foreign_scope_token = session_token.clone();
    foreign_scope_token.ledger_scope = "foreign-scope".to_string();
    assert!(matches!(
        governor.recover_open(&ledger, open_id, foreign_scope_token, None),
        Err(SessionError::MutationConflict {
            kind: "session-open",
            ..
        })
    ));
    let recovered_open = governor
        .recover_open(&ledger, open_id, session_token, None)
        .expect("recover causal open");
    assert_eq!(recovered_open.content_hash(), open_receipt.content_hash());

    assert!(matches!(
        governor.submit_once_durable(&ledger, request_a, "program-a", || {
            panic!("out-of-order recovery invoked caller work")
        }),
        Err(SessionError::RecoveryCausalGap {
            session: 38,
            expected: 1,
            found: 2,
        })
    ));
    assert_eq!(
        governor.consumption(session).unwrap(),
        meter_snapshot_tuple(meter_b.before())
    );
    let replay_b = governor
        .submit_once_durable(&ledger, request_b, "program-b", || {
            panic!("durable request B replay invoked caller work")
        })
        .expect("recover earlier meter commit first");
    assert!(matches!(
        replay_b,
        SubmitOutcome::Duplicate {
            admission_ordinal,
            receipt,
            ref meter_receipt,
            ..
        } if admission_ordinal == admission_b
            && receipt == receipt_b
            && meter_receipt == &meter_b
    ));
    assert_eq!(
        governor.consumption(session).unwrap(),
        meter_snapshot_tuple(meter_b.after())
    );
    let replay_a = governor
        .submit_once_durable(&ledger, request_a, "program-a", || {
            panic!("durable request A replay invoked caller work")
        })
        .expect("recover later meter commit second");
    assert!(matches!(
        replay_a,
        SubmitOutcome::Duplicate {
            admission_ordinal,
            receipt,
            ref meter_receipt,
            ..
        } if admission_ordinal == admission_a
            && receipt == receipt_a
            && meter_receipt == &meter_a
    ));
    assert_eq!(
        governor.consumption(session).unwrap(),
        meter_snapshot_tuple(meter_a.after())
    );
    let no_op = governor
        .flush_scope_to_ledger(&recovered_open.flush_permit(), &ledger)
        .expect("replayed causal state is clean");
    assert_eq!(no_op.committed_terminals, 0);
    assert_eq!(no_op.appended_rows, 0);
    assert!(!no_op.remaining_dirty);
    assert_eq!(
        (
            ledger.table_count("session_claims").unwrap(),
            ledger.table_count("session_terminals").unwrap(),
            ledger.table_count("session_terminal_events").unwrap(),
            ledger.table_count("session_flush_batches").unwrap(),
            ledger.table_count("session_flush_batch_members").unwrap(),
            ledger.table_count("events").unwrap(),
        ),
        durable_counts,
        "durable causal replay changes no registry, witness, or audit row"
    );
}

#[test]
fn ss_003p_pause_acknowledgement_waits_for_pending_submission_meter_commit() {
    let governor = Arc::new(SessionGovernor::new());
    let session = SessionId(39);
    let gate = Arc::new(CancelGate::new());
    let open_id = governor
        .session_open_id(session, "drain-open")
        .expect("open authority");
    governor
        .open_session_gated(open_id, token(39, 100.0, 1e9), Arc::clone(&gate))
        .expect("gated session");
    let request = governor
        .submission_request_id(session, "draining-caller", "draining-program")
        .expect("submission authority");
    let pressure = governor
        .pressure_action_id(session, "draining-pause")
        .expect("pressure authority");
    let refused_while_draining = governor
        .submission_request_id(session, "late-drain-caller", "late-drain-program")
        .expect("same-generation late authority");
    let (started_tx, started_rx) = std::sync::mpsc::channel();
    let (release_tx, release_rx) = std::sync::mpsc::channel();
    let worker = {
        let governor = Arc::clone(&governor);
        std::thread::spawn(move || {
            governor.submit_once(request, || {
                started_tx.send(()).expect("test receiver");
                release_rx.recv().expect("release sender");
                Charge {
                    core_s: 7.0,
                    ..Charge::default()
                }
            })
        })
    };
    started_rx.recv().expect("submission admitted");
    let pressure_receipt = governor
        .apply_memory_pressure(pressure, 3)
        .expect("pause requested");
    let pause_request = pressure_receipt
        .events()
        .iter()
        .find_map(|event| event.pause_request_id)
        .expect("pause request authority");
    let ledger = fs_ledger::Ledger::open(":memory:").expect("fixture ledger");
    let checkpoint = solver_checkpoint(&ledger, pause_request, "settled-checkpoint");
    assert!(matches!(
        governor.submit_once(refused_while_draining, || panic!("draining work ran")),
        Err(SessionError::PauseAlreadyPending { id: 39, .. })
    ));
    assert!(matches!(
        governor.acknowledge_pause(pause_request, &ledger, &checkpoint),
        Err(SessionError::PauseDrainPending {
            id: 39,
            pending_submissions: 1,
        })
    ));

    release_tx.send(()).expect("release submission");
    let executed = worker
        .join()
        .expect("submission worker")
        .expect("submission terminal outcome");
    let meter_receipt = match executed {
        SubmitOutcome::Executed { meter_receipt, .. } => meter_receipt,
        other => panic!("draining submission must execute, got {other:?}"),
    };
    assert_eq!(meter_receipt.after().core_s, 7.0);
    let acknowledgement = governor
        .acknowledge_pause(pause_request, &ledger, &checkpoint)
        .expect("settled generation can rotate");
    assert!(matches!(
        governor
            .submit_once(request, || panic!("terminal replay reran work"))
            .expect("terminal replay"),
        SubmitOutcome::Duplicate {
            meter_receipt: replayed,
            ..
        } if replayed == meter_receipt
    ));
    let refused_while_ready = governor
        .submission_request_id(session, "ready-caller", "ready-program")
        .expect("fresh-generation authority");
    assert!(matches!(
        governor.submit_once(refused_while_ready, || panic!("ready work ran")),
        Err(SessionError::ResumeNotActivated {
            id: 39,
            generation: 1,
        })
    ));
    governor
        .activate_resume(&acknowledgement)
        .expect("fresh generation activated");
    assert!(matches!(
        governor
            .submit_once(refused_while_ready, || Charge {
                core_s: 1.0,
                ..Charge::default()
            })
            .expect("activated generation admits work"),
        SubmitOutcome::Executed { .. }
    ));
    assert_eq!(governor.consumption(session).expect("consumption").0, 8.0);
}

#[test]
#[allow(clippy::too_many_lines)] // One canonical session and scope authority refusal story.
fn ss_003e_ledger_scope_authority_is_canonical_and_fail_closed() {
    let gov = Governor::new();
    for (id, invalid_scope) in [
        (40, ""),
        (41, "main scope"),
        (42, "main\nscope"),
        (43, "máin"),
    ] {
        assert!(matches!(
            gov.open_session(token_in_scope(id, invalid_scope)),
            Err(SessionError::InvalidLedgerScope { .. })
        ));
        assert!(matches!(
            gov.token(SessionId(id)),
            Err(SessionError::UnknownSession { id: unknown }) if unknown == id
        ));
    }
    let oversized = "a".repeat(MAX_LEDGER_SCOPE_BYTES + 1);
    match gov.open_session(token_in_scope(44, &oversized)) {
        Err(SessionError::InvalidLedgerScope {
            scope_preview,
            scope_bytes,
            ..
        }) => {
            assert_eq!(scope_preview.len(), MAX_LEDGER_SCOPE_BYTES);
            assert_eq!(scope_bytes, MAX_LEDGER_SCOPE_BYTES + 1);
        }
        other => panic!("expected bounded oversized-scope refusal, got {other:?}"),
    }
    let split_boundary = format!("{}é", "a".repeat(MAX_LEDGER_SCOPE_BYTES - 1));
    match gov.open_session(token_in_scope(47, &split_boundary)) {
        Err(SessionError::InvalidLedgerScope {
            scope_preview,
            scope_bytes,
            ..
        }) => {
            assert_eq!(scope_preview.len(), MAX_LEDGER_SCOPE_BYTES - 1);
            assert_eq!(scope_bytes, MAX_LEDGER_SCOPE_BYTES + 1);
        }
        other => panic!("expected UTF-8-safe scope preview, got {other:?}"),
    }
    let boundary = "a".repeat(MAX_LEDGER_SCOPE_BYTES);
    let boundary_permit = gov
        .open_session(token_in_scope(44, &boundary))
        .expect("the exact scope-length boundary is admitted after refusal");

    let mut invalid_ops = token_in_scope(49, "ops-invalid");
    invalid_ops.ops = vec!["bad operator".to_string()];
    assert!(matches!(
        gov.open_session(invalid_ops),
        Err(SessionError::InvalidOperatorGrant { index: 0, .. })
    ));
    for (session, malformed) in [(54, "*"), (55, "flux*"), (56, ".flux"), (57, "flux..solve")] {
        let mut malformed_token = token_in_scope(session, "ops-malformed-wildcard");
        malformed_token.ops = vec![malformed.to_string()];
        assert!(matches!(
            gov.open_session(malformed_token),
            Err(SessionError::InvalidOperatorGrant { index: 0, .. })
        ));
    }
    let mut duplicate_ops = token_in_scope(50, "ops-duplicate");
    duplicate_ops.ops = vec!["flux.*".to_string(), "flux.*".to_string()];
    assert!(matches!(
        gov.open_session(duplicate_ops),
        Err(SessionError::DuplicateOperatorGrant { .. })
    ));
    let mut too_many_ops = token_in_scope(51, "ops-count");
    too_many_ops.ops = (0..=MAX_CAPABILITY_OPS)
        .map(|index| format!("op-{index}"))
        .collect();
    assert!(matches!(
        gov.open_session(too_many_ops),
        Err(SessionError::LimitExceeded {
            resource: "capability_operator_grants",
            ..
        })
    ));
    let mut long_op = token_in_scope(52, "ops-bytes");
    long_op.ops = vec!["x".repeat(MAX_CAPABILITY_OP_BYTES + 1)];
    assert!(matches!(
        gov.open_session(long_op),
        Err(SessionError::InvalidOperatorGrant { .. })
    ));
    let mut too_many_op_bytes = token_in_scope(58, "ops-total-bytes");
    too_many_op_bytes.ops = (0..68)
        .map(|index| format!("namespace{index:03}.{}", "x".repeat(110)))
        .collect();
    assert!(matches!(
        gov.open_session(too_many_op_bytes),
        Err(SessionError::LimitExceeded {
            resource: "capability_operator_bytes",
            ..
        })
    ));
    let requested_gate = Arc::new(CancelGate::new());
    requested_gate.request();
    assert!(matches!(
        gov.open_session_gated(
            token_in_scope(53, "pre-requested-gate"),
            Arc::clone(&requested_gate),
        ),
        Err(SessionError::PreRequestedGate { id: 53 })
    ));

    // Reuse an id rejected above: invalid scope admission must not reserve or
    // partially initialize the session.
    let main_permit = gov
        .open_session(token_in_scope(40, "main"))
        .expect("valid scope after atomic refusal");
    gov.charge(
        SessionId(40),
        Charge {
            core_s: 3.0,
            ..Charge::default()
        },
    )
    .expect("dirty main meter");
    let ledger = fs_ledger::Ledger::open(":memory:").expect("scope ledger");
    let foreign_governor = Governor::new();
    let foreign_permit = foreign_governor
        .open_session(token_in_scope(48, "main"))
        .expect("foreign permit fixture");
    assert!(matches!(
        gov.flush_scope_to_ledger(&foreign_permit, &ledger),
        Err(SessionError::ScopePermitMismatch { scope }) if scope == "main"
    ));
    assert!(matches!(
        gov.events_page(&foreign_permit, 0, 1),
        Err(SessionError::ScopePermitMismatch { scope }) if scope == "main"
    ));
    assert!(
        gov.events_page(&main_permit, 0, MAX_EVENT_PAGE_ROWS)
            .is_ok()
    );
    assert!(matches!(
        gov.events_page(&main_permit, 0, MAX_EVENT_PAGE_ROWS + 1),
        Err(SessionError::LimitExceeded {
            resource: "event_page_rows",
            limit: MAX_EVENT_PAGE_ROWS,
            observed_at_least,
        }) if observed_at_least == MAX_EVENT_PAGE_ROWS + 1
    ));
    assert_eq!(ledger.table_count("events").unwrap(), 0);
    gov.flush_scope_to_ledger(&main_permit, &ledger)
        .expect("foreign authority refusal leaves the main cursor dirty");
    assert_eq!(ledger.table_count("events").unwrap(), 2);
    let boundary_ledger = fs_ledger::Ledger::open(":memory:").expect("boundary scope ledger");
    gov.flush_scope_to_ledger(&boundary_permit, &boundary_ledger)
        .expect("main flush did not consume the other scope's meter");
    assert_eq!(boundary_ledger.table_count("events").unwrap(), 1);
}

#[test]
#[allow(clippy::too_many_lines)] // Two-scope cursor/sink/transaction state machine.
fn ss_003f_scoped_flush_isolated_interleaved_retryable_and_sink_bound() {
    const ALPHA: &str = r#"alpha/"quoted"\branch"#;
    const BETA: &str = "beta";
    // Durable protocol (e61t6 fork (a)): a durable governor binds EVERY
    // scope's sink to its recovery ledger at construction, so the old
    // two-sink choreography is impossible by design. Scope isolation is
    // a PER-SCOPE-CURSOR property, not a per-sink one: both scopes
    // flush to the recovery ledger through their own permits, each
    // append covering exactly its own scope's events, and the
    // sink-bound property is proven by foreign-ledger mismatch probes.
    // Submission receipt/failure coverage lives in ss_003c/d.
    let alpha_path = durable_ledger_path("scoped-flush-alpha");
    let alpha_ledger = fs_ledger::Ledger::open(&alpha_path).expect("alpha ledger");
    let gov = Governor::new_durable(&alpha_ledger, DurableGovernorNonce::from_bytes([0xF5; 32]));
    let alpha_permit = gov
        .open_session(token_in_scope(45, ALPHA))
        .expect("canonical JSON-hostile alpha scope");
    let beta_permit = gov
        .open_session(token_in_scope(46, BETA))
        .expect("canonical beta scope");
    gov.flush_scope_to_ledger(&alpha_permit, &alpha_ledger)
        .expect("alpha open prerequisite terminal");
    let alpha_after_open = alpha_ledger.table_count("events").unwrap();
    assert!(matches!(
        gov.submit_once_durable(&alpha_ledger, SessionId(45), "alpha-once", || Charge {
            core_s: 2.0,
            mem_peak_bytes: 5,
            wall_s: 1.0,
        })
        .expect("alpha submission"),
        SubmitOutcome::Executed { .. }
    ));
    gov.apply_memory_pressure(SessionId(45), 1)
        .expect("alpha event one");
    gov.apply_memory_pressure(SessionId(46), 1)
        .expect("interleaved beta event one");
    gov.apply_memory_pressure(SessionId(45), 1)
        .expect("alpha event two");

    let foreign_ledger = fs_ledger::Ledger::open(":memory:").expect("foreign probe ledger");
    alpha_ledger.begin().expect("caller transaction");
    assert!(matches!(
        gov.flush_scope_to_ledger(&alpha_permit, &alpha_ledger),
        Err(SessionError::Persistence { .. })
    ));
    assert_eq!(
        alpha_ledger.table_count("events").unwrap(),
        alpha_after_open
    );
    alpha_ledger.rollback().expect("caller rollback");

    gov.flush_scope_to_ledger(&alpha_permit, &alpha_ledger)
        .expect("alpha retry writes only alpha state");
    let after_alpha_batch = alpha_ledger.table_count("events").unwrap();
    assert_eq!(
        after_alpha_batch,
        alpha_after_open + 3,
        "alpha terminal receipt + two degradation events, nothing of beta"
    );
    gov.flush_scope_to_ledger(&beta_permit, &alpha_ledger)
        .expect("beta flushes through its own cursor into the shared sink");
    let after_beta_batch = alpha_ledger.table_count("events").unwrap();
    assert_eq!(
        after_beta_batch,
        after_alpha_batch + 2,
        "beta open + one degradation event, nothing of alpha re-flushed"
    );
    assert!(
        alpha_ledger.lint().unwrap().is_clean(),
        "the exact quote/backslash-bearing scope must be JSON escaped"
    );

    gov.charge(
        SessionId(45),
        Charge {
            core_s: 1.0,
            ..Charge::default()
        },
    )
    .expect("new alpha meter");
    gov.apply_memory_pressure(SessionId(46), 1)
        .expect("beta event two");
    gov.apply_memory_pressure(SessionId(45), 1)
        .expect("alpha event three");

    assert!(matches!(
        gov.flush_scope_to_ledger(&alpha_permit, &foreign_ledger),
        Err(SessionError::LedgerScopeSinkMismatch { scope, .. }) if scope == ALPHA
    ));
    assert_eq!(
        foreign_ledger.table_count("events").unwrap(),
        0,
        "foreign-sink attempt must append nothing"
    );
    gov.flush_scope_to_ledger(&alpha_permit, &alpha_ledger)
        .expect("wrong-sink attempt leaves both alpha cursors dirty");
    let after_alpha_second = alpha_ledger.table_count("events").unwrap();
    assert_eq!(
        after_alpha_second,
        after_beta_batch + 2,
        "alpha meter + third degradation append after the wrong-sink refusal"
    );
    gov.flush_scope_to_ledger(&beta_permit, &alpha_ledger)
        .expect("alpha activity did not consume beta's degradation cursor");
    let after_beta_second = alpha_ledger.table_count("events").unwrap();
    assert_eq!(after_beta_second, after_alpha_second + 1);

    gov.flush_scope_to_ledger(&alpha_permit, &alpha_ledger)
        .expect("alpha unchanged no-op");
    gov.flush_scope_to_ledger(&beta_permit, &alpha_ledger)
        .expect("beta unchanged no-op");
    assert_eq!(
        alpha_ledger.table_count("events").unwrap(),
        after_beta_second
    );
}

#[test]
fn ss_003g_sink_binding_uses_move_stable_ledger_identity() {
    let gov = Governor::new();
    let permit = gov
        .open_session(token(49, 1e9, 1e9))
        .expect("identity fixture session");
    let ledger = fs_ledger::Ledger::open(":memory:").expect("identity fixture ledger");
    let identity = ledger.instance_id();
    let first = gov
        .flush_scope_to_ledger(&permit, &ledger)
        .expect("bind initial sink");
    assert_eq!(first.appended_rows, 1);

    let ledger = Box::new(ledger);
    assert_eq!(ledger.instance_id(), identity, "moving the handle is inert");
    gov.charge(
        SessionId(49),
        Charge {
            core_s: 1.0,
            ..Charge::default()
        },
    )
    .expect("new dirty meter");
    gov.flush_scope_to_ledger(&permit, &ledger)
        .expect("the moved handle remains the same authority");

    let distinct = fs_ledger::Ledger::open(":memory:").expect("distinct ledger");
    assert_ne!(distinct.instance_id(), identity);
    assert!(matches!(
        gov.flush_scope_to_ledger(&permit, &distinct),
        Err(SessionError::LedgerScopeSinkMismatch {
            bound_sink,
            attempted_sink,
            ..
        }) if bound_sink == identity && attempted_sink == distinct.instance_id()
    ));
}

#[test]
fn ss_003h_bounded_flush_drains_multiple_atomic_chunks_without_duplicates() {
    // Durable protocol (e61t6 fork (a)): the open terminal flushes as
    // a prerequisite, so a trailing degradation event restores the
    // one-row overflow that proves the multi-chunk drain semantics.
    let path = durable_ledger_path("chunk-drain");
    let ledger = fs_ledger::Ledger::open(&path).expect("chunk fixture ledger");
    let gov = Governor::new_durable(&ledger, DurableGovernorNonce::from_bytes([0xA7; 32]));
    let permit = gov
        .open_session(token(50, 1e9, 1e9))
        .expect("chunk fixture session");
    gov.flush_scope_to_ledger(&permit, &ledger)
        .expect("open prerequisite terminal");
    for index in 0..MAX_FLUSH_ROWS {
        assert!(matches!(
            gov.submit_once_durable(
                &ledger,
                SessionId(50),
                &format!("chunk-key-{index:04}"),
                Charge::default
            )
            .expect("bounded key fixture"),
            SubmitOutcome::Executed { .. }
        ));
    }
    gov.apply_memory_pressure(SessionId(50), 1)
        .expect("overflow-row degradation event");
    let first = gov
        .flush_scope_to_ledger(&permit, &ledger)
        .expect("first bounded chunk");
    assert_eq!(first.appended_rows, MAX_FLUSH_ROWS);
    assert!(first.encoded_bytes <= MAX_FLUSH_ENCODED_BYTES);
    assert!(first.remaining_dirty, "one overflow row remains dirty");

    let second = gov
        .flush_scope_to_ledger(&permit, &ledger)
        .expect("second bounded chunk");
    assert_eq!(second.appended_rows, 1);
    assert!(!second.remaining_dirty);
    assert_eq!(
        ledger.table_count("events").unwrap(),
        u64::try_from(MAX_FLUSH_ROWS + 2).expect("fixture count fits"),
        "open + every submission terminal + the degradation event"
    );

    let no_op = gov
        .flush_scope_to_ledger(&permit, &ledger)
        .expect("fully drained flush is a no-op");
    assert_eq!(no_op.appended_rows, 0);
    assert!(!no_op.remaining_dirty);
    assert_eq!(
        ledger.table_count("events").unwrap(),
        u64::try_from(MAX_FLUSH_ROWS + 2).expect("fixture count fits"),
        "cursor commit prevents duplicate rows"
    );
}

#[test]
fn ss_004_estimate_dry_run_and_ledgered_calibration() {
    let node = fs_ir::sexpr::parse(SPOUT).expect("parses");
    let mut models = BTreeMap::new();
    models.insert("xform.level-set-velocity".to_string(), lbm_cost_model());
    let est = estimate(&node, &models, 16.0).expect("valid dry-run inputs");
    assert!(
        (est.wall_p50_s - 409.6).abs() / 409.6 < 0.05,
        "p50 tracks the model: {}",
        est.wall_p50_s
    );
    assert!(est.wall_p10_s <= est.wall_p50_s && est.wall_p50_s <= est.wall_p90_s);
    assert!(est.energy_j > 0.0, "energy estimate present");
    assert_eq!(
        est.mem_ask_bytes,
        Some(96 * 1024 * 1024 * 1024),
        "declared mem ask surfaced"
    );
    assert!(
        est.unmodeled_ops.contains(&"ascent.optimize".to_string()),
        "coverage gaps are stated, not silent"
    );
    let undotted = fs_ir::sexpr::parse("(study \"undotted-model\" (simulate :size 4))")
        .expect("undotted operation parses");
    models.insert("simulate".to_string(), lbm_cost_model());
    let undotted_est = estimate(&undotted, &models, 1.0)
        .expect("a registered undotted operation is modeled by registry identity");
    assert!(undotted_est.wall_p50_s > 0.0);
    assert!(
        undotted_est.unmodeled_ops.is_empty(),
        "a registered undotted operation must not disappear from dry-run coverage"
    );
    // Calibration: synthetic actuals at 1.1x the estimate.
    let calib = CalibrationReport::new();
    for k in 0..20 {
        let mut e = est.clone();
        e.wall_p50_s *= 1.0 + f64::from(k) * 0.01;
        calib
            .record(&e, e.wall_p50_s * 1.1)
            .expect("finite calibration row");
    }
    let (q10, q50, q90) = calib.ratio_quantiles().expect("rows");
    assert!((q50 - 1.1).abs() < 1e-9, "median ratio is the true bias");
    assert!(q10 <= q50 && q50 <= q90);
    // Ledgered as a content-addressed artifact.
    let dir = std::env::temp_dir().join(format!("fs-session-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("mkdir");
    let ledger =
        fs_ledger::Ledger::open(dir.join("cal.led").to_str().expect("utf8")).expect("ledger");
    let hash = calib.flush_to_ledger(&ledger).expect("flush");
    let bytes = ledger.get_artifact(&hash).expect("get").expect("present");
    let text = String::from_utf8(bytes).expect("utf8");
    assert!(text.contains("estimate-calibration") && text.contains("ratio_quantiles"));
    let _ = std::fs::remove_dir_all(&dir);
    verdict(
        "ss-004",
        "dry-run p10/p50/p90 + energy + honest coverage; calibration ledgered",
    );
}

#[test]
#[allow(clippy::too_many_lines)] // One hostile resource-domain conformance matrix.
fn ss_004b_estimate_refuses_invalid_resource_domains() {
    let valid = fs_ir::sexpr::parse(SPOUT).expect("valid fixture");
    let mut models = BTreeMap::new();
    models.insert("xform.level-set-velocity".to_string(), lbm_cost_model());
    for cores in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1.0] {
        assert!(matches!(
            estimate(&valid, &models, cores),
            Err(SessionError::InvalidResource {
                resource: "estimate cores",
                ..
            })
        ));
    }

    for memory in [
        "-1GiB",
        "0B",
        "0.1B",
        "0.99999999999999999B",
        "1.00000000000000001B",
        "1e-1B",
        "100000000000000000000GiB",
        "3cores",
    ] {
        let study = fs_ir::sexpr::parse(&format!(
            "(study \"invalid-estimate\" (budget (mem {memory})))"
        ))
        .expect("the IR can represent the hostile Count fixture");
        assert!(matches!(
            estimate(&study, &models, 1.0),
            Err(SessionError::InvalidResource {
                resource: "declared memory ask",
                ..
            })
        ));
    }
    for (memory, expected) in [
        ("1.5GiB", 3_u64 << 29),
        ("1e3B", 1_000),
        ("18446744073709551615B", u64::MAX),
    ] {
        let study = fs_ir::sexpr::parse(&format!(
            "(study \"exact-estimate\" (budget (mem {memory})))"
        ))
        .expect("exact memory fixture parses");
        assert_eq!(
            estimate(&study, &models, 1.0)
                .expect("exact memory estimate")
                .mem_ask_bytes,
            Some(expected),
            "{memory} changed authority during estimation"
        );
    }
    for malformed in [
        "(study \"missing-memory\" (budget (mem)))",
        "(study \"bare-memory\" (budget mem))",
        "(study \"keyword-memory\" (budget :mem 1GiB))",
        "(study \"headless-memory\" (budget (:mem 1GiB)))",
        "(study \"empty-budget-entry\" (budget ()))",
        "(study \"wrong-memory-kind\" (budget (mem 3kg)))",
        "(study \"extra-memory-operand\" (budget (mem 1GiB 2GiB)))",
        "(study \"duplicate-memory\" (budget (mem 1GiB) (mem 2GiB)))",
    ] {
        let study = fs_ir::sexpr::parse(malformed).expect("malformed budget shape remains an AST");
        assert!(matches!(
            estimate(&study, &models, 1.0),
            Err(SessionError::Submission { .. })
        ));
    }
    let body_decoy = fs_ir::sexpr::parse("(study \"body-decoy\" (budget (wall 1s)) (mem 7GiB))")
        .expect("body decoy parses");
    assert_eq!(
        estimate(&body_decoy, &models, 1.0)
            .expect("a body call named mem is not a budget declaration")
            .mem_ask_bytes,
        None
    );
    let budget_before_decoy =
        fs_ir::sexpr::parse("(study \"budget-before-decoy\" (budget (mem 1GiB)) (mem 7GiB))")
            .expect("budget and body decoy parse");
    assert_eq!(
        estimate(&budget_before_decoy, &models, 1.0)
            .expect("only the recognized budget supplies memory")
            .mem_ask_bytes,
        Some(1024 * 1024 * 1024)
    );
    let negative_size =
        fs_ir::sexpr::parse("(study \"invalid-size\" (xform.level-set-velocity field :dof -1))")
            .expect("negative operation sizes are representable before admission");
    assert!(matches!(
        estimate(&negative_size, &models, 1.0),
        Err(SessionError::InvalidResource {
            resource: "estimate operation size",
            ..
        })
    ));
    for malformed_size in [
        "(study \"missing-size\" (xform.level-set-velocity field :dof))",
        "(study \"nonnumeric-size\" (xform.level-set-velocity field :size \"many\"))",
        "(study \"duplicate-size\" (xform.level-set-velocity field :dof 4 :size 8))",
    ] {
        let study = fs_ir::sexpr::parse(malformed_size).expect("malformed size remains an AST");
        assert!(matches!(
            estimate(&study, &models, 1.0),
            Err(SessionError::Submission { .. })
        ));
    }
    let implicit_unit_size =
        fs_ir::sexpr::parse("(study \"implicit-unit-size\" (xform.level-set-velocity field))")
            .expect("size-free operation remains a valid AST");
    assert!(
        estimate(&implicit_unit_size, &models, 1.0).is_ok(),
        "the unit-size default remains available only when no size feature is declared"
    );

    let calibration = CalibrationReport::new();
    let finite = Estimate {
        wall_p10_s: 1.0,
        wall_p50_s: 1.0,
        wall_p90_s: 1.0,
        mem_ask_bytes: None,
        energy_j: 45.0,
        unmodeled_ops: Vec::new(),
        weakest_cost_evidence: None,
    };
    for actual in [f64::NAN, f64::INFINITY, -1.0] {
        assert!(calibration.record(&finite, actual).is_err());
    }
    let ratio_overflow = Estimate {
        wall_p50_s: f64::MIN_POSITIVE,
        ..finite.clone()
    };
    assert!(calibration.record(&ratio_overflow, f64::MAX).is_err());
    assert!(calibration.ratio_quantiles().is_none());
    assert_eq!(
        calibration.to_json(),
        "{\"kind\":\"estimate-calibration\",\"rows\":[],\"ratio_quantiles\":null,\
         \"zero_predictions\":{\"true_zero\":0,\"unmodeled\":0,\"actual_quantiles_s\":null}}"
    );
}

#[test]
fn ss_005_degradation_ladder_declared_order_and_pause_resume() {
    #[derive(Debug, PartialEq)]
    struct ToySolver {
        step: u64,
        field: Vec<f64>,
    }
    impl SolverState for ToySolver {
        const TYPE_ID: u64 = 0x544f_5953_4f4c_0001;
        const SCHEMA_VERSION: u32 = 1;

        fn encode(&self, enc: &mut codec::Enc) {
            enc.put_u64(self.step);
            enc.put_u64(self.field.len() as u64);
            for v in &self.field {
                enc.put_f64(*v);
            }
        }
        fn decode(dec: &mut codec::Dec<'_>) -> Result<Self, codec::CodecError> {
            let step = dec.get_u64()?;
            let n = usize::try_from(dec.get_u64()?).expect("fits");
            let mut field = Vec::with_capacity(n);
            for _ in 0..n {
                field.push(dec.get_f64()?);
            }
            Ok(ToySolver { step, field })
        }
    }
    let gov = Governor::new();
    let gate = Arc::new(CancelGate::new());
    let permit = gov
        .open_session_gated(token(5, 1e9, 1e9), Arc::clone(&gate))
        .expect("valid token");
    // Level 1: only the spill step fires.
    let l1 = gov.apply_memory_pressure(SessionId(5), 1).expect("session");
    assert_eq!(l1.len(), 1);
    assert_eq!(l1[0].step, DegradationStep::SpillColdArenas);
    assert!(!gate.is_requested(), "level 1 must not pause");
    // Level 3: all three fire IN THE DECLARED ORDER; pause requests the
    // session's OWN gate (bound at open — no gate crosses the API).
    let l3 = gov.apply_memory_pressure(SessionId(5), 3).expect("session");
    let steps: Vec<DegradationStep> = l3.iter().map(|e| e.step).collect();
    assert_eq!(
        steps,
        vec![
            DegradationStep::SpillColdArenas,
            DegradationStep::CoarsenAdaptively,
            DegradationStep::PauseSerializeResume
        ],
        "the ladder order is the contract"
    );
    assert!(gate.is_requested(), "pause step requests cancellation");
    // Pause-serialize-resume equality (P7): snapshot round-trips exactly.
    let solver = ToySolver {
        step: 4242,
        field: (0..64).map(|i| f64::from(i) * 0.25 - 3.0).collect(),
    };
    let bytes = solver.to_bytes();
    let resumed = ToySolver::from_bytes(&bytes).expect("resume");
    assert_eq!(resumed, solver, "pause-serialize-resume must be lossless");
    // Events are attributed and ordinal-ordered.
    let events = gov
        .events_page(&permit, 0, MAX_EVENT_PAGE_ROWS)
        .expect("bounded event page");
    assert_eq!(events.len(), 4);
    assert!(events.windows(2).all(|w| w[0].ordinal < w[1].ordinal));
    assert!(events.iter().all(|e| !e.attribution.is_empty()));
    verdict(
        "ss-005",
        "ladder fires spill->coarsen->pause in declared order; snapshot round-trip exact",
    );
}

#[test]
#[allow(clippy::too_many_lines)] // One owned-gate pause lifecycle conformance scenario.
fn ss_011_pressure_actions_bind_to_owned_session_gates() {
    // Bead gp3.13 acceptance battery: gates are OWNED (bound at open),
    // wrong-session pauses unrepresentable, out-of-ladder levels fail,
    // and a pause is never complete without a generation-bound ledger receipt.
    let gov = Governor::new();
    let ledger = fs_ledger::Ledger::open(":memory:").expect("fixture ledger");
    let gate_a = Arc::new(CancelGate::new());
    let gate_b = Arc::new(CancelGate::new());
    let permit_a = gov
        .open_session_gated(token(41, 1e9, 1e9), Arc::clone(&gate_a))
        .expect("valid token");
    gov.open_session_gated(token(42, 1e9, 1e9), Arc::clone(&gate_b))
        .expect("valid token");
    gov.open_session(token(43, 1e9, 1e9)).expect("valid token"); // ungated

    // (a) Levels 0 and > 3 are REFUSED, never clamped, nothing ledgered.
    for bad in [0u8, 4, 200] {
        assert_eq!(
            gov.apply_memory_pressure(SessionId(41), bad),
            Err(SessionError::InvalidPressureLevel { level: bad }),
            "level {bad} must be refused"
        );
    }
    assert!(
        gov.events_page(&permit_a, 0, MAX_EVENT_PAGE_ROWS)
            .expect("bounded event page")
            .is_empty(),
        "refused levels must not ledger events"
    );

    // (b) Level 3 on an UNGATED session is refused ATOMICALLY: no gate
    // to reach the computation means no pause claim and no partial
    // ladder application.
    assert_eq!(
        gov.apply_memory_pressure(SessionId(43), 3),
        Err(SessionError::UngatedSession { id: 43 }),
        "ungated session must refuse level 3"
    );
    assert!(
        gov.events_page(&permit_a, 0, MAX_EVENT_PAGE_ROWS)
            .expect("bounded event page")
            .is_empty(),
        "a refused pause must not half-apply the ladder"
    );
    // Levels 1-2 need no gate: spill/coarsen are synchronous.
    let l2 = gov
        .apply_memory_pressure(SessionId(43), 2)
        .expect("levels 1-2 need no gate");
    assert_eq!(l2.len(), 2);
    assert!(l2.iter().all(|e| e.phase == StepPhase::Declared));

    // (c) Level 3 on session A requests ONLY A's gate; B is untouched.
    let l3 = gov.apply_memory_pressure(SessionId(41), 3).expect("gated");
    assert!(gate_a.is_requested(), "the target session's gate fires");
    assert!(
        !gate_b.is_requested(),
        "level 3 must request only the target session"
    );
    // (d) The request event is phase Requested — NOT complete.
    let pause = l3.last().expect("three steps fired");
    assert_eq!(pause.step, DegradationStep::PauseSerializeResume);
    assert_eq!(pause.phase, StepPhase::Requested);
    let request_id = pause
        .pause_request_id
        .expect("level three mints opaque request authority");
    assert_eq!(request_id.session(), SessionId(41));
    assert_eq!(request_id.gate_generation(), 0);
    assert_eq!(request_id.requested_ordinal(), pause.ordinal);
    assert!(gov.pause_pending(SessionId(41)).expect("known session"));

    let before_repeat = gov
        .events_page(&permit_a, 0, MAX_EVENT_PAGE_ROWS)
        .expect("bounded event page");
    assert_eq!(
        gov.apply_memory_pressure(SessionId(41), 3),
        Err(SessionError::PauseAlreadyPending {
            id: 41,
            requested_ordinal: pause.ordinal,
        }),
        "a repeated level-3 request must refuse before replaying spill/coarsen"
    );
    assert_eq!(
        gov.apply_memory_pressure(SessionId(41), 1),
        Err(SessionError::PauseAlreadyPending {
            id: 41,
            requested_ordinal: pause.ordinal,
        }),
        "synchronous pressure also refuses while a checkpoint generation is draining"
    );
    assert_eq!(
        gov.events_page(&permit_a, 0, MAX_EVENT_PAGE_ROWS)
            .expect("bounded event page"),
        before_repeat,
        "repeated level-3 refusal must not mutate the scoped event stream"
    );

    // (e) A live old worker prevents the executor from minting drain proof.
    let authority = request_id.checkpoint_authority();
    let mut run_bytes = [0_u8; 8];
    run_bytes.copy_from_slice(&authority.as_bytes()[..8]);
    let run = RunId(u64::from_le_bytes(run_bytes));
    let drain_gate = CancelGate::new_clock_free();
    let tracker = DrainTracker::new(run, &drain_gate);
    let old_worker = tracker.register_worker().expect("old worker admitted");
    drain_gate.request();
    assert!(matches!(
        tracker.finalize(),
        Err(fs_exec::DrainFinalizeError::WorkersStillRunning { active: 1 })
    ));
    assert!(gov.pause_pending(SessionId(41)).expect("known session"));
    drop(old_worker);
    let report = tracker.finalize().expect("old worker drained");
    let snapshot =
        fs_exec::solver::envelope::seal(0x4653_434b_5054, 1, run.0, b"solver-state-0xf00d");
    let artifact = ledger
        .put_artifact(
            fs_ledger::session_registry::SOLVER_STATE_ARTIFACT_KIND,
            &snapshot,
            None,
        )
        .expect("solver-state artifact");
    let checkpoint_claim = fs_ledger::session_registry::SolverCheckpointClaim {
        session: request_id.session().0,
        pause_authority: authority,
        gate_generation: request_id.gate_generation(),
        solver_state_artifact: artifact.hash,
        drain_report: &report,
    };
    let checkpoint = ledger
        .attest_solver_checkpoint(checkpoint_claim)
        .expect("checkpoint receipt");
    assert_eq!(
        ledger
            .attest_solver_checkpoint(checkpoint_claim)
            .expect("response-loss retry"),
        checkpoint,
        "retry after losing the mint response must recover one receipt"
    );
    let recovered_checkpoint = ledger
        .solver_checkpoint_receipt(authority)
        .expect("checkpoint lookup")
        .expect("stored checkpoint");
    assert_eq!(recovered_checkpoint, checkpoint);

    // Malformed transport cannot forge a typed receipt.
    let mut forged_transport = checkpoint.to_bytes();
    forged_transport[20] ^= 1;
    assert!(
        fs_ledger::session_registry::SolverCheckpointReceipt::from_bytes(&forged_transport)
            .is_err()
    );

    // (f) Cross-session, foreign-ledger, and cross-governor authorities fail closed.
    let foreign = Governor::new();
    foreign
        .open_session_gated(token(42, 1e9, 1e9), Arc::new(CancelGate::new()))
        .expect("foreign fixture session");
    let foreign_request = foreign
        .apply_memory_pressure(SessionId(42), 3)
        .expect("foreign pause request")
        .last()
        .and_then(|event| event.pause_request_id)
        .expect("foreign request authority");
    let cross_session_ledger = fs_ledger::Ledger::open(":memory:").expect("foreign ledger");
    let cross_session_checkpoint =
        solver_checkpoint(&cross_session_ledger, foreign_request, "cross-session");
    assert!(matches!(
        gov.acknowledge_pause(request_id, &cross_session_ledger, &cross_session_checkpoint,),
        Err(SessionError::PauseCheckpointMismatch {
            id: 41,
            reason: "cross-session",
            ..
        })
    ));
    assert!(matches!(
        gov.acknowledge_pause(
            foreign_request,
            &cross_session_ledger,
            &cross_session_checkpoint,
        ),
        Err(SessionError::PauseRequestMismatch { id: 42, .. })
    ));
    let foreign_ledger = fs_ledger::Ledger::open(":memory:").expect("foreign ledger");
    let foreign_checkpoint = solver_checkpoint(&foreign_ledger, request_id, "foreign-ledger");
    assert!(matches!(
        gov.acknowledge_pause(request_id, &ledger, &foreign_checkpoint),
        Err(SessionError::PauseCheckpointMismatch {
            id: 41,
            reason: "unverified-ledger-receipt",
            ..
        })
    ));

    // (g) The verified receipt is the ONLY route to Complete; the
    // completion event cites the request it acknowledges.
    let acknowledged = gov
        .acknowledge_pause(request_id, &ledger, &recovered_checkpoint)
        .expect("pending pause");
    let done = acknowledged.event();
    assert_eq!(done.phase, StepPhase::Complete);
    assert_eq!(done.gate_generation, Some(0));
    assert_eq!(acknowledged.resume_generation(), 1);
    assert!(!acknowledged.resume_gate().is_requested());
    let checkpoint = done
        .checkpoint
        .as_ref()
        .expect("complete event carries structured checkpoint evidence");
    assert_eq!(
        checkpoint.preview(),
        recovered_checkpoint.content_hash().to_hex()
    );
    assert_eq!(checkpoint.byte_len(), 64);
    assert!(
        done.attribution
            .contains(&recovered_checkpoint.content_hash().to_string())
    );
    assert!(done.attribution.contains(&artifact.hash.to_string()));
    assert!(
        done.attribution
            .contains(&report.content_hash().to_string())
    );
    assert!(
        done.attribution
            .contains(&format!("ordinal {}", pause.ordinal)),
        "completion must cite the request it acknowledges"
    );
    assert!(!gov.pause_pending(SessionId(41)).expect("known session"));
    // Lost-response replay is idempotent and recovers the exact same gate.
    let replayed = gov
        .acknowledge_pause(request_id, &ledger, &recovered_checkpoint)
        .expect("identical acknowledgement replay");
    assert_eq!(replayed.event(), acknowledged.event());
    assert_eq!(
        replayed.resume_generation(),
        acknowledged.resume_generation()
    );
    assert!(Arc::ptr_eq(
        &replayed.resume_gate(),
        &acknowledged.resume_gate()
    ));
    // Foreign-ledger evidence cannot replace the completed generation.
    assert!(matches!(
        gov.acknowledge_pause(request_id, &ledger, &foreign_checkpoint),
        Err(SessionError::PauseCheckpointMismatch { id: 41, .. })
    ));

    // No pressure transition can request the fresh gate before resumed workers
    // explicitly adopt it. If an external owner requests that never-activated
    // gate, activation refuses and identical acknowledgement replay replaces
    // only the Arc identity at the SAME generation.
    assert!(matches!(
        gov.apply_memory_pressure(SessionId(41), 3),
        Err(SessionError::ResumeNotActivated {
            id: 41,
            generation: 1,
        })
    ));
    let events_before_gate_recovery = gov
        .events_page(&permit_a, 0, MAX_EVENT_PAGE_ROWS)
        .expect("bounded pre-recovery event page");
    let cancelled_before_activation = acknowledged.resume_gate();
    cancelled_before_activation.request();
    assert_eq!(
        gov.activate_resume(&acknowledged),
        Err(SessionError::ResumeGateAlreadyRequested {
            id: 41,
            generation: 1,
        })
    );
    let recovered = gov
        .acknowledge_pause(request_id, &ledger, &recovered_checkpoint)
        .expect("identical replay replaces a never-activated requested gate");
    assert_eq!(recovered.event(), acknowledged.event());
    assert_eq!(recovered.resume_generation(), 1);
    assert!(!recovered.resume_gate().is_requested());
    assert!(!Arc::ptr_eq(
        &recovered.resume_gate(),
        &cancelled_before_activation
    ));
    assert_eq!(
        gov.activate_resume(&acknowledged),
        Err(SessionError::ResumeAcknowledgementMismatch { id: 41 }),
        "the acknowledgement carrying the replaced Arc must stay stale"
    );
    let recovered_replay = gov
        .acknowledge_pause(request_id, &ledger, &recovered_checkpoint)
        .expect("replacement acknowledgement replay");
    assert!(Arc::ptr_eq(
        &recovered.resume_gate(),
        &recovered_replay.resume_gate()
    ));
    assert_eq!(
        gov.events_page(&permit_a, 0, MAX_EVENT_PAGE_ROWS)
            .expect("bounded post-recovery event page"),
        events_before_gate_recovery,
        "replacing a never-activated gate is not another pause generation"
    );
    assert!(matches!(
        gov.apply_memory_pressure(SessionId(41), 1),
        Err(SessionError::ResumeNotActivated {
            id: 41,
            generation: 1,
        })
    ));
    gov.activate_resume(&recovered)
        .expect("workers adopted fresh gate");
    gov.activate_resume(&recovered)
        .expect("activation replay is idempotent");

    // A second cycle requests the fresh generation rather than inheriting the
    // monotonic cancellation state of the drained first gate.
    let second = gov
        .apply_memory_pressure(SessionId(41), 3)
        .expect("fresh gate admits a second pause cycle");
    let second_pause = second.last().expect("second pause request");
    assert_eq!(second_pause.gate_generation, Some(1));
    let second_request = second_pause
        .pause_request_id
        .expect("second request authority");
    assert!(recovered.resume_gate().is_requested());
    gov.activate_resume(&recovered)
        .expect("activation replay remains idempotent after the next pause request");
    assert!(gov.pause_pending(SessionId(41)).expect("known session"));
    assert!(matches!(
        gov.acknowledge_pause(second_request, &ledger, &recovered_checkpoint),
        Err(SessionError::PauseCheckpointMismatch {
            id: 41,
            reason: "stale-generation",
            ..
        })
    ));
    let second_checkpoint = solver_checkpoint(&ledger, second_request, "second-generation");
    let second_ack = gov
        .acknowledge_pause(second_request, &ledger, &second_checkpoint)
        .expect("second pending pause");
    assert_eq!(second_ack.event().gate_generation, Some(1));
    assert_eq!(second_ack.resume_generation(), 2);
    assert!(!second_ack.resume_gate().is_requested());
    assert!(matches!(
        gov.acknowledge_pause(request_id, &ledger, &recovered_checkpoint),
        Err(SessionError::PauseRequestMismatch { id: 41, .. })
    ));
    // (h) The ledgered stream never contains an unacknowledged Complete:
    // exactly one Complete, and it follows its Requested ordinal.
    let events = gov
        .events_page(&permit_a, 0, MAX_EVENT_PAGE_ROWS)
        .expect("bounded event page");
    let completes: Vec<_> = events
        .iter()
        .filter(|e| e.phase == StepPhase::Complete)
        .collect();
    assert_eq!(completes.len(), 2);
    assert!(completes[0].ordinal > pause.ordinal);
    assert!(events.windows(2).all(|w| w[0].ordinal < w[1].ordinal));
    verdict(
        "ss-011",
        "pressure binds to owned gates: bad levels refused, ungated level-3 atomic refusal, \
         target-only request, complete only via a drained ledger checkpoint receipt",
    );
}

#[test]
fn ss_006_budget_infeasible_surfaces_as_ranked_guidance() {
    // The §11.3 canonical fixture: admission's BudgetInfeasible finding
    // becomes a Guidance value with cost-model-ranked fixes.
    let src = SPOUT.replace("(wall 2h)", "(wall 60s)");
    let node = fs_ir::sexpr::parse(&src).expect("parses");
    let mut cost_models = BTreeMap::new();
    cost_models.insert("xform.level-set-velocity".to_string(), lbm_cost_model());
    let cx = fs_ir::admission::AdmissionContext {
        cost_freshness: None,
        router: None,
        chart_requirements: Vec::new(),
        cost_models,
        capability: Some(
            token(9, 1e9, 1e9)
                .to_admission()
                .expect("bounded grants project"),
        ),
        regime: None,
        regime_policy: fs_ir::admission::RegimePolicy::Warn,
    };
    let report = fs_ir::admission::admit(&node, &cx);
    assert!(!report.admitted);
    // The sealed cost-model doctrine (2pmb) adds a provisional-evidence
    // WARN under the same "budget" check before the rejection, so the
    // guidance fixture must select the REJECT finding specifically.
    let finding = report
        .findings
        .iter()
        .find(|f| f.check == "budget" && f.severity == fs_ir::admission::Severity::Reject)
        .expect("budget rejection finding");
    let guidance = Guidance::from_finding(finding);
    assert_eq!(guidance.code, "budget-rejection");
    assert!(guidance.diagnosis.contains("BudgetInfeasible"));
    assert!(
        guidance.fixes.len() >= 2,
        "ranked fixes travel with the refusal"
    );
    let rendered = guidance.render();
    assert!(rendered.contains("fix#0") && rendered.contains("predicted wall"));
    verdict(
        "ss-006",
        "BudgetInfeasible teaches: code + diagnosis + ranked fixes render",
    );
}

#[test]
fn ss_007_governor_storm_structured_outcomes_only() {
    let gov = Arc::new(Governor::new());
    for id in 0..8u64 {
        // Adversarial: tiny grants on odd sessions.
        let grant = if id % 2 == 0 { 1e6 } else { 20.0 };
        gov.open_session(token(id, grant, 1e9))
            .expect("valid token");
    }
    let mut handles = Vec::new();
    for id in 0..8u64 {
        for worker in 0..4u32 {
            let gov = Arc::clone(&gov);
            handles.push(std::thread::spawn(move || {
                let key = format!("storm:{id}:{worker}");
                let out = gov
                    .submit_once(SessionId(id), &key, || Charge {
                        core_s: 9.0,
                        ..Charge::default()
                    })
                    .expect("known session");
                let enforce = gov
                    .charge(
                        SessionId(id),
                        Charge {
                            core_s: 9.0,
                            ..Charge::default()
                        },
                    )
                    .expect("known session");
                (out, enforce)
            }));
        }
    }
    let mut throttled_or_paused = 0usize;
    for h in handles {
        let (out, enforce) = h.join().expect("join");
        assert!(
            matches!(out, SubmitOutcome::Executed { .. }),
            "unique keys execute"
        );
        match enforce {
            Enforcement::Ok => {}
            Enforcement::Throttled { .. } | Enforcement::Paused { .. } => {
                throttled_or_paused += 1;
            }
        }
    }
    assert!(
        throttled_or_paused > 0,
        "adversarial grants must trip enforcement somewhere"
    );
    // Every session's meters are exact: 4 submits + 4 charges x 9 core-s.
    for id in 0..8u64 {
        let (core_s, ..) = gov.consumption(SessionId(id)).expect("meters");
        assert!((core_s - 72.0).abs() < 1e-9, "session {id}: {core_s}");
    }
    verdict(
        "ss-007",
        "32-way storm over adversarial grants: exact meters, structured outcomes only",
    );
}

#[test]
fn ss_008_panicking_submission_is_non_blocking_and_terminal() {
    const WAITERS: usize = 8;
    let gov = Arc::new(Governor::new());
    gov.open_session(token(80, 1e9, 1e9)).expect("valid token");
    let executions = Arc::new(AtomicU32::new(0));
    let (started_tx, started_rx) = std::sync::mpsc::channel();
    let (release_tx, release_rx) = std::sync::mpsc::channel();
    let key = "panic-terminal".to_string();
    let panic_text = "seeded submission panic".to_string();

    let owner = {
        let gov = Arc::clone(&gov);
        let executions = Arc::clone(&executions);
        let key = key.clone();
        let panic_text = panic_text.clone();
        std::thread::spawn(move || {
            gov.submit_once(SessionId(80), &key, || {
                executions.fetch_add(1, Ordering::SeqCst);
                started_tx.send(()).expect("test receiver alive");
                release_rx.recv().expect("release sender alive");
                std::panic::panic_any(panic_text);
            })
        })
    };

    started_rx
        .recv_timeout(std::time::Duration::from_secs(5))
        .expect("owner reached caller work after installing Pending");
    let mut waiters = Vec::new();
    for _ in 0..WAITERS {
        let gov = Arc::clone(&gov);
        let key = key.clone();
        waiters.push(std::thread::spawn(move || {
            gov.submit_once(SessionId(80), &key, || {
                panic!("a duplicate must never execute");
            })
            .expect("session remains valid")
        }));
    }
    for waiter in waiters {
        assert_eq!(
            waiter.join().expect("waiter returned"),
            SubmitOutcome::InFlight,
            "a duplicate of Pending must return immediately without blocking"
        );
    }
    release_tx.send(()).expect("owner is still waiting");
    let owner_outcome = owner
        .join()
        .expect("panic was contained by submit_once")
        .expect("session remains valid");

    assert_eq!(executions.load(Ordering::SeqCst), 1, "work ran once");
    let (admission_ordinal, receipt, evidence) = match owner_outcome {
        SubmitOutcome::Failed {
            admission_ordinal,
            receipt,
            evidence,
        } => (admission_ordinal, receipt, evidence),
        other => panic!("owner must publish one terminal failure, got {other:?}"),
    };
    assert_eq!(evidence.preview(), panic_text);
    assert_eq!(evidence.byte_len(), panic_text.len());
    let request_id = gov
        .inner
        .submission_request_id(SessionId(80), &key, &key)
        .expect("typed panic request");
    assert!(receipt.matches_failure(request_id, "main", admission_ordinal, &evidence));
    assert_eq!(
        gov.consumption(SessionId(80)).expect("meters").0.to_bits(),
        0.0f64.to_bits(),
        "failed work is never charged"
    );

    let retry = gov
        .submit_once(SessionId(80), &key, || {
            executions.fetch_add(1, Ordering::SeqCst);
            Charge::default()
        })
        .expect("terminal failure is readable");
    assert!(
        matches!(retry, SubmitOutcome::Failed { receipt: r, evidence: ref e, .. } if r == receipt && e == &evidence),
        "same-key retry receives the terminal failure; a new key is explicit retry"
    );
    assert_eq!(executions.load(Ordering::SeqCst), 1);
    assert!(matches!(
        gov.submit_once(SessionId(80), "panic-terminal:retry-1", Charge::default)
            .expect("new key executes"),
        SubmitOutcome::Executed { .. }
    ));
    verdict(
        "ss-008",
        "seeded panic: one execution, bounded full-digest terminal evidence, Pending callers return InFlight, no charge",
    );
}

#[test]
fn ss_009_invalid_resources_fail_closed_without_poisoning_meters() {
    let gov = Governor::new();
    let invalid = [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1.0];
    let mut next_id = 90u64;
    for value in invalid {
        for field in ["core_s", "wall_s"] {
            let id = next_id;
            next_id += 1;
            let mut bad = token(id, 100.0, 100.0);
            match field {
                "core_s" => bad.core_s = value,
                "wall_s" => bad.wall_s = value,
                _ => unreachable!(),
            }
            assert!(matches!(
                gov.open_session(bad),
                Err(SessionError::InvalidResource { .. })
            ));
            assert!(
                matches!(
                    gov.token(SessionId(id)),
                    Err(SessionError::UnknownSession { .. })
                ),
                "invalid {field} grant must not partially open session {id}"
            );
        }
    }

    gov.open_session(token(102, 100.0, 1_000.0))
        .expect("valid token");
    let before = gov.consumption(SessionId(102)).expect("meters");
    for value in invalid {
        for delta in [
            Charge {
                core_s: value,
                ..Charge::default()
            },
            Charge {
                wall_s: value,
                ..Charge::default()
            },
        ] {
            assert!(matches!(
                gov.charge(SessionId(102), delta),
                Err(SessionError::InvalidResource { .. })
            ));
            assert_eq!(
                gov.consumption(SessionId(102)).expect("meters"),
                before,
                "invalid charge must not mutate any meter"
            );
        }
    }

    verdict(
        "ss-009",
        "NaN/infinite/negative grants and deltas rejected before mutation",
    );
}

#[test]
fn ss_010_exact_grant_throttles_and_accumulated_overflow_is_atomic() {
    let gov = Governor::new();
    gov.open_session(token(102, 100.0, 1_000.0))
        .expect("valid token");
    let at_grant = gov
        .charge(
            SessionId(102),
            Charge {
                core_s: 100.0,
                ..Charge::default()
            },
        )
        .expect("valid charge");
    assert!(
        matches!(
            at_grant,
            Enforcement::Throttled { used, granted, .. }
                if used.to_bits() == granted.to_bits()
        ),
        "the documented exact-grant boundary throttles"
    );
    assert!(matches!(
        gov.charge(
            SessionId(102),
            Charge {
                core_s: 20.0,
                ..Charge::default()
            }
        )
        .expect("valid charge"),
        Enforcement::Throttled { .. }
    ));
    assert!(matches!(
        gov.charge(
            SessionId(102),
            Charge {
                core_s: 1.0,
                ..Charge::default()
            }
        )
        .expect("valid charge"),
        Enforcement::Paused { .. }
    ));

    gov.open_session(token(103, f64::MAX, f64::MAX))
        .expect("finite maximum grants are valid");
    let _ = gov
        .charge(
            SessionId(103),
            Charge {
                core_s: f64::MAX,
                ..Charge::default()
            },
        )
        .expect("first finite charge");
    assert!(matches!(
        gov.charge(
            SessionId(103),
            Charge {
                core_s: f64::MAX,
                ..Charge::default()
            }
        ),
        Err(SessionError::InvalidResource { .. })
    ));
    assert_eq!(
        gov.consumption(SessionId(103)).expect("meters").0.to_bits(),
        f64::MAX.to_bits(),
        "overflow rejection leaves the prior finite meter intact"
    );
    verdict(
        "ss-010",
        "exact grant throttles; accumulated overflow is refused without mutating the finite meter",
    );
}

#[test]
fn ss_012_session_and_idempotency_collection_caps_are_exact_and_atomic() {
    let sessions = Governor::new();
    for id in 0..MAX_SESSIONS_PER_SCOPE {
        sessions
            .open_session(token_in_scope(
                u64::try_from(id).expect("fixture id fits"),
                "crowded",
            ))
            .expect("exact per-scope boundary is admitted");
    }
    let refused_id = u64::try_from(MAX_SESSIONS_PER_SCOPE).expect("fixture id fits");
    assert!(matches!(
        sessions.open_session(token_in_scope(refused_id, "crowded")),
        Err(SessionError::LimitExceeded {
            resource: "sessions_per_scope",
            limit: MAX_SESSIONS_PER_SCOPE,
            observed_at_least,
        }) if observed_at_least == MAX_SESSIONS_PER_SCOPE + 1
    ));
    assert!(matches!(
        sessions.token(SessionId(refused_id)),
        Err(SessionError::UnknownSession { id }) if id == refused_id
    ));

    for id in MAX_SESSIONS_PER_SCOPE..MAX_SESSIONS_PER_GOVERNOR {
        let id = u64::try_from(id).expect("fixture id fits");
        sessions
            .open_session(token_in_scope(id, &format!("scope-{id}")))
            .expect("fill exact governor boundary");
    }
    let governor_overflow = u64::try_from(MAX_SESSIONS_PER_GOVERNOR).expect("fixture id fits");
    assert!(matches!(
        sessions.open_session(token_in_scope(governor_overflow, "overflow")),
        Err(SessionError::LimitExceeded {
            resource: "sessions_per_governor",
            limit: MAX_SESSIONS_PER_GOVERNOR,
            observed_at_least,
        }) if observed_at_least == MAX_SESSIONS_PER_GOVERNOR + 1
    ));

    let keys = Governor::new();
    keys.open_session(token(9_000, 1e9, 1e9))
        .expect("key-cap fixture session");
    for index in 0..MAX_IDEMPOTENCY_KEYS_PER_SESSION {
        assert!(matches!(
            keys.submit_once(
                SessionId(9_000),
                &format!("bounded-key-{index:04}"),
                Charge::default,
            )
            .expect("exact key boundary is admitted"),
            SubmitOutcome::Executed { .. }
        ));
    }
    let overflow_work = AtomicU32::new(0);
    assert!(matches!(
        keys.submit_once(SessionId(9_000), "one-key-too-many", || {
            overflow_work.fetch_add(1, Ordering::SeqCst);
            Charge::default()
        }),
        Err(SessionError::LimitExceeded {
            resource: "idempotency_keys_per_session",
            limit: MAX_IDEMPOTENCY_KEYS_PER_SESSION,
            observed_at_least,
        }) if observed_at_least == MAX_IDEMPOTENCY_KEYS_PER_SESSION + 1
    ));
    assert_eq!(overflow_work.load(Ordering::SeqCst), 0);
    assert!(matches!(
        keys.submit_once(SessionId(9_000), "bounded-key-0000", || {
            panic!("terminal duplicate at the cap must not rerun")
        })
        .expect("existing terminal keys remain queryable"),
        SubmitOutcome::Duplicate { .. }
    ));
}

#[test]
#[allow(clippy::too_many_lines)] // One real-reopen proof for the complete typed lifecycle chain.
fn ss_013_durable_meter_and_l3_lifecycle_reopen_without_state_or_row_drift() {
    let path = durable_ledger_path("lifecycle");
    let nonce = DurableGovernorNonce::from_bytes([0xD3; 32]);
    let ledger = fs_ledger::Ledger::open(&path).expect("on-disk lifecycle ledger");
    let governor = SessionGovernor::new_durable(&ledger, nonce).expect("durable governor");
    let session = SessionId(9_013);
    let token = token_in_scope(session.0, "durable-lifecycle");
    let initial_gate = Arc::new(CancelGate::new());
    let open_id = governor
        .session_open_id(session, "durable-lifecycle-open")
        .expect("open authority");
    let open_receipt = governor
        .open_session_gated(open_id, token.clone(), Arc::clone(&initial_gate))
        .expect("gated open");
    let companion_session = SessionId(9_014);
    let companion_token = token_in_scope(companion_session.0, "durable-lifecycle");
    let companion_open_id = governor
        .session_open_id(companion_session, "durable-companion-open")
        .expect("companion open authority");
    governor
        .open_session(companion_open_id, companion_token.clone())
        .expect("companion open");
    let permit = open_receipt.flush_permit();
    governor
        .flush_scope_to_ledger(&permit, &ledger)
        .expect("open prerequisite terminal");

    let delta = Charge {
        core_s: 17.0,
        mem_peak_bytes: 19,
        wall_s: 23.0,
    };
    let meter_id = governor
        .meter_report_id(session, "durable-meter")
        .expect("meter authority");
    let meter_receipt = governor.charge(meter_id, delta).expect("meter commit");
    let l1_id = governor
        .pressure_action_id(session, "durable-pressure-l1")
        .expect("L1 authority");
    let l1 = governor.apply_memory_pressure(l1_id, 1).expect("L1 action");
    let l2_id = governor
        .pressure_action_id(session, "durable-pressure-l2")
        .expect("L2 authority");
    let l2 = governor.apply_memory_pressure(l2_id, 2).expect("L2 action");
    let l3_id = governor
        .pressure_action_id(session, "durable-pressure-l3")
        .expect("L3 authority");
    let l3 = governor.apply_memory_pressure(l3_id, 3).expect("L3 action");
    assert!(initial_gate.is_requested());
    let pause_request = l3
        .events()
        .iter()
        .find_map(|event| event.pause_request_id)
        .expect("L3 request authority");
    let companion_action_id = governor
        .pressure_action_id(companion_session, "interleaved-companion-l1")
        .expect("interleaved companion action authority");
    let companion_l1 = governor
        .apply_memory_pressure(companion_action_id, 1)
        .expect("interleaved companion action");
    let checkpoint = solver_checkpoint(&ledger, pause_request, "durable-checkpoint");
    let acknowledgement = governor
        .acknowledge_pause(pause_request, &ledger, &checkpoint)
        .expect("pause acknowledgement");
    let activation = governor
        .activate_resume(&acknowledgement)
        .expect("resume activation");
    let resumed_submission_id = governor
        .submission_request_id(session, "resumed-slot", "resumed-program")
        .expect("resumed-generation submission authority");
    let resumed_submission = governor
        .submit_once_durable(&ledger, resumed_submission_id, "resumed-program", || {
            Charge {
                core_s: 29.0,
                ..Charge::default()
            }
        })
        .expect("durable submission after gate rotation");
    let resumed_submission_receipt = match resumed_submission {
        SubmitOutcome::Executed { receipt, .. } => receipt,
        other => panic!("expected resumed execution, got {other:?}"),
    };
    let second_l3_id = governor
        .pressure_action_id(session, "durable-pressure-l3-second")
        .expect("second L3 authority");
    let second_l3 = governor
        .apply_memory_pressure(second_l3_id, 3)
        .expect("second L3 request after activation");
    let second_pause_request = second_l3
        .events()
        .iter()
        .find_map(|event| event.pause_request_id)
        .expect("second L3 request authority");
    let second_checkpoint =
        solver_checkpoint(&ledger, second_pause_request, "durable-checkpoint-second");
    let second_acknowledgement = governor
        .acknowledge_pause(second_pause_request, &ledger, &second_checkpoint)
        .expect("second pause acknowledgement");
    let flushed = governor
        .flush_scope_to_ledger(&permit, &ledger)
        .expect("meter and lifecycle terminal batch");
    assert_eq!(flushed.committed_terminals, 10);
    assert_eq!(flushed.appended_rows, 14);
    assert!(!flushed.remaining_dirty);
    let counts = (
        ledger.table_count("session_claims").unwrap(),
        ledger.table_count("session_terminals").unwrap(),
        ledger.table_count("session_terminal_events").unwrap(),
        ledger.table_count("session_flush_batches").unwrap(),
        ledger.table_count("session_flush_batch_members").unwrap(),
        ledger.table_count("events").unwrap(),
    );
    assert_eq!(counts.0, 12);
    assert_eq!(counts.1, 12);
    assert_eq!(counts.2, 16);
    assert_eq!(counts.5, 16);
    drop(governor);
    drop(ledger);

    let ledger = fs_ledger::Ledger::open(&path).expect("reopened lifecycle ledger");
    let governor = SessionGovernor::new_durable(&ledger, nonce).expect("reopened governor");
    let recovered_initial_gate = Arc::new(CancelGate::new());
    let recovered_open = governor
        .recover_open(
            &ledger,
            open_id,
            token,
            Some(Arc::clone(&recovered_initial_gate)),
        )
        .expect("recover gated open");
    governor
        .recover_open(&ledger, companion_open_id, companion_token.clone(), None)
        .expect("recover companion open");
    assert_eq!(recovered_open.content_hash(), open_receipt.content_hash());
    assert_eq!(
        governor
            .recover_meter(&ledger, meter_id, delta)
            .expect("recover meter"),
        meter_receipt
    );
    assert_eq!(
        governor
            .recover_pressure(&ledger, l1_id, 1)
            .expect("recover L1"),
        l1
    );
    assert_eq!(
        governor
            .recover_pressure(&ledger, l2_id, 2)
            .expect("recover L2"),
        l2
    );
    assert_eq!(
        governor
            .recover_pressure(&ledger, l3_id, 3)
            .expect("recover L3"),
        l3
    );
    assert!(recovered_initial_gate.is_requested());
    let recovered_resume_gate = Arc::new(CancelGate::new());
    assert!(matches!(
        governor.recover_pause_acknowledgement(
            &ledger,
            pause_request,
            &checkpoint,
            Arc::clone(&recovered_resume_gate),
        ),
        Err(SessionError::TerminalCorrupt { .. })
    ));
    assert_eq!(
        governor
            .recover_pressure(&ledger, companion_action_id, 1)
            .expect("recover interleaved companion action"),
        companion_l1
    );
    let recovered_acknowledgement = governor
        .recover_pause_acknowledgement(
            &ledger,
            pause_request,
            &checkpoint,
            Arc::clone(&recovered_resume_gate),
        )
        .expect("recover pause acknowledgement");
    assert_eq!(
        recovered_acknowledgement.content_hash(),
        acknowledgement.content_hash()
    );
    assert_eq!(recovered_acknowledgement.event(), acknowledgement.event());
    assert!(Arc::ptr_eq(
        &recovered_acknowledgement.resume_gate(),
        &recovered_resume_gate
    ));
    assert_eq!(
        governor
            .recover_resume_activation(&ledger, &recovered_acknowledgement)
            .expect("recover activation"),
        activation
    );
    let executions = AtomicU32::new(0);
    assert!(matches!(
        governor
            .submit_once_durable(
                &ledger,
                resumed_submission_id,
                "resumed-program",
                || {
                    executions.fetch_add(1, Ordering::SeqCst);
                    Charge::default()
                },
            )
            .expect("recover resumed-generation submission"),
        SubmitOutcome::Duplicate { receipt, .. } if receipt == resumed_submission_receipt
    ));
    assert_eq!(executions.load(Ordering::SeqCst), 0);
    assert_eq!(
        governor
            .recover_pressure(&ledger, second_l3_id, 3)
            .expect("recover next L3 request"),
        second_l3
    );
    assert_eq!(
        governor
            .recover_resume_activation(&ledger, &recovered_acknowledgement)
            .expect("prior activation remains replayable after next gate request"),
        activation
    );

    assert_eq!(
        governor
            .recover_pressure(&ledger, l3_id, 3)
            .expect("terminal pressure replay after activation"),
        l3
    );
    assert_eq!(
        governor
            .recover_pause_acknowledgement(
                &ledger,
                pause_request,
                &checkpoint,
                Arc::clone(&recovered_resume_gate),
            )
            .expect("exact acknowledgement replay")
            .content_hash(),
        acknowledgement.content_hash()
    );
    assert_eq!(
        governor
            .recover_resume_activation(&ledger, &recovered_acknowledgement)
            .expect("exact activation replay"),
        activation
    );
    let recovered_second_gate = Arc::new(CancelGate::new());
    assert_eq!(
        governor
            .recover_pause_acknowledgement(
                &ledger,
                second_pause_request,
                &second_checkpoint,
                Arc::clone(&recovered_second_gate),
            )
            .expect("recover second acknowledgement")
            .content_hash(),
        second_acknowledgement.content_hash()
    );
    assert!(matches!(
        governor.recover_pause_acknowledgement(
            &ledger,
            pause_request,
            &checkpoint,
            Arc::clone(&recovered_resume_gate),
        ),
        Err(SessionError::PauseAcknowledgementConflict { .. })
    ));
    assert!(matches!(
        governor.recover_pause_acknowledgement(
            &ledger,
            pause_request,
            &checkpoint,
            Arc::new(CancelGate::new()),
        ),
        Err(SessionError::PauseAcknowledgementConflict { .. })
    ));
    let foreign_checkpoint_ledger =
        fs_ledger::Ledger::open(":memory:").expect("foreign checkpoint ledger");
    let altered_checkpoint = solver_checkpoint(
        &foreign_checkpoint_ledger,
        pause_request,
        "altered-checkpoint",
    );
    assert!(matches!(
        governor.recover_pause_acknowledgement(
            &ledger,
            pause_request,
            &altered_checkpoint,
            Arc::clone(&recovered_resume_gate),
        ),
        Err(SessionError::PauseCheckpointMismatch { .. })
    ));
    let foreign = fs_ledger::Ledger::open(":memory:").expect("foreign ledger");
    assert!(matches!(
        governor.recover_meter(&foreign, meter_id, delta),
        Err(SessionError::RecoveryLedgerMismatch { .. })
    ));
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

    // Recovery order is intentionally lifecycle-first here. Authenticated
    // terminal meter/submission receipts from older generations must remain
    // installable later in their own contiguous meter order.
    let ledger = fs_ledger::Ledger::open(&path).expect("historical-order ledger");
    let governor = SessionGovernor::new_durable(&ledger, nonce).expect("historical governor");
    let historical_initial_gate = Arc::new(CancelGate::new());
    governor
        .recover_open(
            &ledger,
            open_id,
            token_in_scope(session.0, "durable-lifecycle"),
            Some(Arc::clone(&historical_initial_gate)),
        )
        .expect("historical-order open");
    governor
        .recover_open(&ledger, companion_open_id, companion_token, None)
        .expect("historical companion open");
    governor
        .recover_pressure(&ledger, l1_id, 1)
        .expect("historical L1");
    governor
        .recover_pressure(&ledger, l2_id, 2)
        .expect("historical L2");
    governor
        .recover_pressure(&ledger, l3_id, 3)
        .expect("historical first L3");
    governor
        .recover_pressure(&ledger, companion_action_id, 1)
        .expect("historical interleaved companion action");
    let historical_gate_one = Arc::new(CancelGate::new());
    let historical_ack_one = governor
        .recover_pause_acknowledgement(
            &ledger,
            pause_request,
            &checkpoint,
            Arc::clone(&historical_gate_one),
        )
        .expect("historical first acknowledgement");
    governor
        .recover_resume_activation(&ledger, &historical_ack_one)
        .expect("historical first activation");
    governor
        .recover_pressure(&ledger, second_l3_id, 3)
        .expect("historical second L3");
    let historical_gate_two = Arc::new(CancelGate::new());
    governor
        .recover_pause_acknowledgement(
            &ledger,
            second_pause_request,
            &second_checkpoint,
            historical_gate_two,
        )
        .expect("historical second acknowledgement");
    assert_eq!(
        governor
            .recover_meter(&ledger, meter_id, delta)
            .expect("late historical meter"),
        meter_receipt
    );
    let historical_executions = AtomicU32::new(0);
    assert!(matches!(
        governor
            .submit_once_durable(
                &ledger,
                resumed_submission_id,
                "resumed-program",
                || {
                    historical_executions.fetch_add(1, Ordering::SeqCst);
                    Charge::default()
                },
            )
            .expect("late historical submission"),
        SubmitOutcome::Duplicate { receipt, .. } if receipt == resumed_submission_receipt
    ));
    assert_eq!(historical_executions.load(Ordering::SeqCst), 0);
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
    verdict(
        "ss-013",
        "durable meter and complete L3 lifecycle recover exactly after real ledger reopen",
    );
}

#[test]
fn ss_004a_estimate_carries_weakest_cost_evidence() {
    // Sealed authority (bead 2pmb): weakest-wins, never upgraded.
    let node = fs_ir::sexpr::parse(SPOUT).expect("parses");
    let empty: BTreeMap<String, SealedCostModel> = BTreeMap::new();
    let unmodeled_est = estimate(&node, &empty, 4.0).expect("estimates without models");
    assert_eq!(
        unmodeled_est.weakest_cost_evidence, None,
        "no modeled calls means no evidence claim at all"
    );
    let mut models = BTreeMap::new();
    models.insert("xform.level-set-velocity".to_string(), lbm_cost_model());
    let est = estimate(&node, &models, 4.0).expect("estimates");
    assert_eq!(
        est.weakest_cost_evidence,
        Some(fs_plan::CostEvidenceClass::ProvisionalUnaudited),
        "one provisional contributor marks the whole estimate"
    );
}

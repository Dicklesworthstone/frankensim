//! fs-session conformance (the gp3.7 bead). Acceptance: budget
//! enforcement throttles/pauses with structured outcomes (never a silent
//! kill); double-submit with one idempotency key = one execution, one
//! charge (concurrency-stress-tested); estimate() accuracy tracked vs
//! actuals with a ledgered calibration report; the degradation ladder
//! fires in its declared order under synthetic memory pressure with
//! pause-serialize-resume equality; errors surface as ranked guidance.

use fs_exec::CancelGate;
use fs_exec::solver::{SolverState, codec};
use fs_plan::{CostModel, CostObservation};
use fs_session::{
    CalibrationReport, CapabilityToken, Charge, DegradationStep, Enforcement, Estimate, Governor,
    Guidance, MAX_EVENT_PAGE_ROWS, MAX_FLUSH_ENCODED_BYTES, MAX_FLUSH_ROWS,
    MAX_IDEMPOTENCY_KEYS_PER_SESSION, MAX_LEDGER_SCOPE_BYTES, MAX_RETAINED_EVIDENCE_BYTES,
    MAX_SESSIONS_PER_GOVERNOR, MAX_SESSIONS_PER_SCOPE, SessionError, SessionId, StepPhase,
    SubmitOutcome, estimate,
};
use std::collections::BTreeMap;
use std::sync::Arc;
use std::sync::atomic::{AtomicU32, Ordering};

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

const SPOUT: &str = r#"(study "spout-laminar-v3"
  (seed 0x5EED0001) (versions (constellation :lock "2026-07"))
  (budget (wall 2h) (mem 96GiB) (qoi-rel-error 2e-2))
  (let lever (xform.level-set-velocity vessel :band 12mm :dof 4096))
  (ascent.optimize J :over lever :method (lbfgs :m 17)))"#;

fn lbm_cost_model() -> CostModel {
    let obs: Vec<CostObservation> = (1..=12)
        .map(|k| CostObservation {
            size: f64::from(k) * 512.0,
            cost_s: 0.1 * f64::from(k) * 512.0,
        })
        .collect();
    CostModel::fit(&obs).expect("fits")
}

#[test]
fn ss_001_token_bridges_into_static_admission() {
    let t = token(1, 3600.0, 7200.0);
    assert!(t.grants_op("flux.free-surface-lbm"));
    assert!(!t.grants_op("quantum.anneal"));
    let cap = t.to_admission();
    assert!((cap.wall_s - 7200.0).abs() < f64::EPSILON);
    assert_eq!(cap.mem_bytes, t.mem_bytes);
    assert_eq!(cap.cores, t.cores);

    let mut boundary = token(2, 1.0, 1.0);
    boundary.mem_bytes = u64::MAX;
    boundary.cores = 9_007_199_254_740_993;
    let boundary_cap = boundary.to_admission();
    assert_eq!(boundary_cap.mem_bytes, u64::MAX);
    assert_eq!(boundary_cap.cores, 9_007_199_254_740_993);
    // The bridge feeds fs-ir admission directly.
    let node = fs_ir::sexpr::parse(SPOUT).expect("parses");
    let cx = fs_ir::admission::AdmissionContext {
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
    let key = Governor::idempotency_key("agent-a", SPOUT);
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
    let key = Governor::idempotency_key("agent:alpha", "program:beta");
    assert!(key.starts_with("fs-session-idem-v2:"));
    assert_eq!(key.len(), "fs-session-idem-v2:".len() + 64);
    assert_eq!(
        key,
        Governor::idempotency_key("agent:alpha", "program:beta")
    );
    assert_ne!(
        key,
        Governor::idempotency_key("agent", "alpha:program:beta"),
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
    let positive = opened(10.0)
        .submit_once(SessionId(33), key, || positive_zero_charge)
        .expect("positive-zero submit");
    let (positive_receipt, positive_enforcement) = match positive {
        SubmitOutcome::Executed {
            receipt,
            enforcement,
            ..
        } => (receipt, enforcement),
        other => panic!("expected execution, got {other:?}"),
    };
    assert_eq!(positive_enforcement, Enforcement::Ok);
    assert!(positive_receipt.matches_success(
        SessionId(33),
        "main",
        key,
        positive_zero_charge,
        &positive_enforcement,
    ));
    assert!(
        !positive_receipt.matches_success(
            SessionId(33),
            "other-scope",
            key,
            positive_zero_charge,
            &positive_enforcement,
        ),
        "receipt v2 binds the immutable ledger scope"
    );
    assert!(
        !positive_receipt.matches_success(
            SessionId(33),
            "main",
            key,
            Charge {
                mem_peak_bytes: positive_zero_charge.mem_peak_bytes + 1,
                ..positive_zero_charge
            },
            &positive_enforcement,
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

    let failed = opened(10.0)
        .submit_once(SessionId(33), key, || panic!("receipt-bound failure"))
        .expect("panic is a terminal outcome");
    let failed_receipt = match failed {
        SubmitOutcome::Failed { receipt, evidence } => {
            assert!(receipt.matches_failure(SessionId(33), "main", key, &evidence));
            assert!(evidence.preview().contains("receipt-bound failure"));
            receipt
        }
        other => panic!("expected failed receipt, got {other:?}"),
    };
    assert_ne!(positive_receipt, failed_receipt);

    let throttled = Governor::new();
    let throttled_permit = throttled
        .open_session(token(33, 1.0, 1e9))
        .expect("valid receipt-test token");
    let charge = Charge {
        core_s: 1.0,
        ..Charge::default()
    };
    let executed = throttled
        .submit_once(SessionId(33), key, || charge)
        .expect("throttled work still completes");
    let (receipt, enforcement) = match executed {
        SubmitOutcome::Executed {
            receipt,
            enforcement,
            ..
        } => (receipt, enforcement),
        other => panic!("expected execution, got {other:?}"),
    };
    assert!(matches!(enforcement, Enforcement::Throttled { .. }));
    assert!(receipt.matches_success(SessionId(33), "main", key, charge, &enforcement));
    assert!(matches!(
        throttled
            .submit_once(SessionId(33), key, || panic!("duplicate ran"))
            .expect("duplicate"),
        SubmitOutcome::Duplicate {
            receipt: duplicate_receipt,
            enforcement: Enforcement::Throttled { .. },
        } if duplicate_receipt == receipt
    ));

    let ledger = fs_ledger::Ledger::open(":memory:").expect("receipt ledger");
    throttled
        .flush_scope_to_ledger(&throttled_permit, &ledger)
        .expect("strict JSON receipt event");
    assert_eq!(ledger.table_count("events").unwrap(), 2);
    assert!(ledger.lint().unwrap().is_clean());
}

#[test]
fn ss_003d_ledger_flush_is_atomic_incremental_and_retryable() {
    let gov = Governor::new();
    let permit = gov
        .open_session(token(34, 1e9, 1e9))
        .expect("flush-test session");
    let key = "flush-once";
    assert!(matches!(
        gov.submit_once(SessionId(34), key, || Charge {
            core_s: 2.0,
            mem_peak_bytes: 5,
            wall_s: 1.0,
        })
        .expect("submission"),
        SubmitOutcome::Executed { .. }
    ));
    let ledger = fs_ledger::Ledger::open(":memory:").expect("flush ledger");

    ledger.begin().expect("caller transaction");
    let refused = gov
        .flush_scope_to_ledger(&permit, &ledger)
        .expect_err("flush cannot promise durability inside a caller transaction");
    assert!(matches!(refused, SessionError::Persistence { .. }));
    assert_eq!(ledger.table_count("events").unwrap(), 0);
    ledger.rollback().expect("caller rollback");

    gov.flush_scope_to_ledger(&permit, &ledger)
        .expect("the refused batch remains fully dirty for retry");
    assert_eq!(
        ledger.table_count("events").unwrap(),
        2,
        "one consumption snapshot and one terminal idempotency event"
    );
    gov.flush_scope_to_ledger(&permit, &ledger)
        .expect("unchanged no-op flush");
    assert_eq!(
        ledger.table_count("events").unwrap(),
        2,
        "repeated flush must not duplicate semantic events"
    );
    assert!(matches!(
        gov.submit_once(SessionId(34), key, || panic!("duplicate ran"))
            .expect("terminal duplicate"),
        SubmitOutcome::Duplicate { .. }
    ));
    gov.flush_scope_to_ledger(&permit, &ledger)
        .expect("duplicate observation is not a new event");
    assert_eq!(ledger.table_count("events").unwrap(), 2);

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
    assert_eq!(ledger.table_count("events").unwrap(), 3);
    gov.apply_memory_pressure(SessionId(34), 1)
        .expect("one degradation event");
    gov.flush_scope_to_ledger(&permit, &ledger)
        .expect("new degradation event appends once");
    assert_eq!(ledger.table_count("events").unwrap(), 4);
    gov.flush_scope_to_ledger(&permit, &ledger)
        .expect("second no-op flush");
    assert_eq!(ledger.table_count("events").unwrap(), 4);
    assert!(ledger.lint().unwrap().is_clean());
}

#[test]
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
    assert_eq!(ledger.table_count("events").unwrap(), 1);
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
    let gov = Governor::new();
    let alpha_permit = gov
        .open_session(token_in_scope(45, ALPHA))
        .expect("canonical JSON-hostile alpha scope");
    let beta_permit = gov
        .open_session(token_in_scope(46, BETA))
        .expect("canonical beta scope");
    assert!(matches!(
        gov.submit_once(SessionId(45), "alpha-once", || Charge {
            core_s: 2.0,
            mem_peak_bytes: 5,
            wall_s: 1.0,
        })
        .expect("alpha submission"),
        SubmitOutcome::Executed { .. }
    ));
    assert!(matches!(
        gov.submit_once(SessionId(46), "beta-once", || Charge {
            core_s: 7.0,
            mem_peak_bytes: 11,
            wall_s: 4.0,
        })
        .expect("beta submission"),
        SubmitOutcome::Executed { .. }
    ));
    assert!(matches!(
        gov.submit_once(SessionId(46), r#"beta-"failed"\key"#, || Charge {
            core_s: f64::NAN,
            ..Charge::default()
        })
        .expect("invalid charge becomes one terminal failure receipt"),
        SubmitOutcome::Failed { .. }
    ));
    gov.apply_memory_pressure(SessionId(45), 1)
        .expect("alpha event one");
    gov.apply_memory_pressure(SessionId(46), 1)
        .expect("interleaved beta event one");
    gov.apply_memory_pressure(SessionId(45), 1)
        .expect("alpha event two");

    let alpha_ledger = fs_ledger::Ledger::open(":memory:").expect("alpha ledger");
    let beta_ledger = fs_ledger::Ledger::open(":memory:").expect("beta ledger");
    alpha_ledger.begin().expect("caller transaction");
    assert!(matches!(
        gov.flush_scope_to_ledger(&alpha_permit, &alpha_ledger),
        Err(SessionError::Persistence { .. })
    ));
    assert_eq!(alpha_ledger.table_count("events").unwrap(), 0);
    alpha_ledger.rollback().expect("caller rollback");

    gov.flush_scope_to_ledger(&alpha_permit, &alpha_ledger)
        .expect("alpha retry writes only alpha state");
    assert_eq!(
        alpha_ledger.table_count("events").unwrap(),
        4,
        "alpha meter + terminal receipt + two alpha degradation events"
    );
    assert_eq!(beta_ledger.table_count("events").unwrap(), 0);
    gov.flush_scope_to_ledger(&beta_permit, &beta_ledger)
        .expect("beta independently binds its own sink");
    assert_eq!(
        beta_ledger.table_count("events").unwrap(),
        4,
        "beta meter + success/failure receipts + one beta degradation event"
    );
    assert!(
        alpha_ledger.lint().unwrap().is_clean(),
        "the exact quote/backslash-bearing scope must be JSON escaped"
    );
    assert!(beta_ledger.lint().unwrap().is_clean());

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
        gov.flush_scope_to_ledger(&alpha_permit, &beta_ledger),
        Err(SessionError::LedgerScopeSinkMismatch { scope, .. }) if scope == ALPHA
    ));
    assert_eq!(
        beta_ledger.table_count("events").unwrap(),
        4,
        "cross-scope sink attempt must append nothing"
    );
    gov.flush_scope_to_ledger(&alpha_permit, &alpha_ledger)
        .expect("wrong-sink attempt leaves both alpha cursors dirty");
    assert_eq!(alpha_ledger.table_count("events").unwrap(), 6);
    gov.flush_scope_to_ledger(&beta_permit, &beta_ledger)
        .expect("alpha activity did not consume beta's degradation cursor");
    assert_eq!(beta_ledger.table_count("events").unwrap(), 5);

    gov.flush_scope_to_ledger(&alpha_permit, &alpha_ledger)
        .expect("alpha unchanged no-op");
    gov.flush_scope_to_ledger(&beta_permit, &beta_ledger)
        .expect("beta unchanged no-op");
    assert_eq!(alpha_ledger.table_count("events").unwrap(), 6);
    assert_eq!(beta_ledger.table_count("events").unwrap(), 5);
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
    let gov = Governor::new();
    let permit = gov
        .open_session(token(50, 1e9, 1e9))
        .expect("chunk fixture session");
    for index in 0..MAX_FLUSH_ROWS {
        assert!(matches!(
            gov.submit_once(
                SessionId(50),
                &format!("chunk-key-{index:04}"),
                Charge::default
            )
            .expect("bounded key fixture"),
            SubmitOutcome::Executed { .. }
        ));
    }
    let ledger = fs_ledger::Ledger::open(":memory:").expect("chunk fixture ledger");
    let first = gov
        .flush_scope_to_ledger(&permit, &ledger)
        .expect("first bounded chunk");
    assert_eq!(first.appended_rows, MAX_FLUSH_ROWS);
    assert!(first.encoded_bytes <= MAX_FLUSH_ENCODED_BYTES);
    assert!(first.remaining_dirty, "one terminal row remains dirty");

    let second = gov
        .flush_scope_to_ledger(&permit, &ledger)
        .expect("second bounded chunk");
    assert_eq!(second.appended_rows, 1);
    assert!(!second.remaining_dirty);
    assert_eq!(
        ledger.table_count("events").unwrap(),
        u64::try_from(MAX_FLUSH_ROWS + 1).expect("fixture count fits")
    );

    let no_op = gov
        .flush_scope_to_ledger(&permit, &ledger)
        .expect("fully drained flush is a no-op");
    assert_eq!(no_op.appended_rows, 0);
    assert!(!no_op.remaining_dirty);
    assert_eq!(
        ledger.table_count("events").unwrap(),
        u64::try_from(MAX_FLUSH_ROWS + 1).expect("fixture count fits"),
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
    // and a pause is never complete without a checkpoint receipt.
    let gov = Governor::new();
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
    assert!(l2.iter().all(|e| e.phase == StepPhase::Applied));

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
        gov.events_page(&permit_a, 0, MAX_EVENT_PAGE_ROWS)
            .expect("bounded event page"),
        before_repeat,
        "repeated level-3 refusal must not mutate the scoped event stream"
    );

    // (e) A blank receipt cannot complete the pause (and the pending
    // request survives the refusal).
    assert!(matches!(
        gov.acknowledge_pause(SessionId(41), "  "),
        Err(SessionError::Submission { .. })
    ));
    assert!(gov.pause_pending(SessionId(41)).expect("known session"));
    // (f) Acknowledging a session with NO outstanding request is refused.
    assert_eq!(
        gov.acknowledge_pause(SessionId(42), "ckpt-b-1"),
        Err(SessionError::NoPendingPause { id: 42 })
    );
    // (g) The checkpoint receipt is the ONLY route to Complete; the
    // completion event cites the request it acknowledges.
    let checkpoint_receipt = format!(
        "solver-state-0xf00d:{}",
        "é".repeat(MAX_RETAINED_EVIDENCE_BYTES)
    );
    let done = gov
        .acknowledge_pause(SessionId(41), &checkpoint_receipt)
        .expect("pending pause");
    assert_eq!(done.phase, StepPhase::Complete);
    assert!(done.attribution.contains("solver-state-0xf00d"));
    let checkpoint = done
        .checkpoint
        .as_ref()
        .expect("complete event carries structured checkpoint evidence");
    assert_eq!(checkpoint.preview().len(), MAX_RETAINED_EVIDENCE_BYTES);
    assert_eq!(checkpoint.byte_len(), checkpoint_receipt.len());
    assert!(
        done.attribution
            .contains(&format!("ordinal {}", pause.ordinal)),
        "completion must cite the request it acknowledges"
    );
    assert!(!gov.pause_pending(SessionId(41)).expect("known session"));
    // Double-acknowledgement is refused (the claim is consumed).
    assert_eq!(
        gov.acknowledge_pause(SessionId(41), &checkpoint_receipt),
        Err(SessionError::NoPendingPause { id: 41 })
    );
    // (h) The ledgered stream never contains an unacknowledged Complete:
    // exactly one Complete, and it follows its Requested ordinal.
    let events = gov
        .events_page(&permit_a, 0, MAX_EVENT_PAGE_ROWS)
        .expect("bounded event page");
    let completes: Vec<_> = events
        .iter()
        .filter(|e| e.phase == StepPhase::Complete)
        .collect();
    assert_eq!(completes.len(), 1);
    assert!(completes[0].ordinal > pause.ordinal);
    assert!(events.windows(2).all(|w| w[0].ordinal < w[1].ordinal));
    verdict(
        "ss-011",
        "pressure binds to owned gates: bad levels refused, ungated level-3 atomic refusal, \
         target-only request, complete only via checkpoint receipt",
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
        router: None,
        chart_requirements: Vec::new(),
        cost_models,
        capability: Some(token(9, 1e9, 1e9).to_admission()),
        regime: None,
        regime_policy: fs_ir::admission::RegimePolicy::Warn,
    };
    let report = fs_ir::admission::admit(&node, &cx);
    assert!(!report.admitted);
    let finding = report
        .findings
        .iter()
        .find(|f| f.check == "budget")
        .expect("budget finding");
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
    let (receipt, evidence) = match owner_outcome {
        SubmitOutcome::Failed { receipt, evidence } => (receipt, evidence),
        other => panic!("owner must publish one terminal failure, got {other:?}"),
    };
    assert_eq!(evidence.preview(), panic_text);
    assert_eq!(evidence.byte_len(), panic_text.len());
    assert!(receipt.matches_failure(SessionId(80), "main", &key, &evidence));
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
        matches!(retry, SubmitOutcome::Failed { receipt: r, evidence: ref e } if r == receipt && e == &evidence),
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

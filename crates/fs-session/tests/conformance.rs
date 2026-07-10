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
    CalibrationReport, CapabilityToken, Charge, DegradationStep, Enforcement, Governor, Guidance,
    SessionError, SessionId, SubmitOutcome, estimate,
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
        cores: 16.0,
        ledger_scope: "main".to_string(),
    }
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
        if let SubmitOutcome::Duplicate { receipt: r } = o {
            assert_eq!(*r, receipt, "duplicates share the original receipt");
        }
    }
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
        "16-thread race: one execution, one charge, shared receipt",
    );
}

#[test]
fn ss_004_estimate_dry_run_and_ledgered_calibration() {
    let node = fs_ir::sexpr::parse(SPOUT).expect("parses");
    let mut models = BTreeMap::new();
    models.insert("xform.level-set-velocity".to_string(), lbm_cost_model());
    let est = estimate(&node, &models, 16.0);
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
    // Calibration: synthetic actuals at 1.1x the estimate.
    let calib = CalibrationReport::new();
    for k in 0..20 {
        let mut e = est.clone();
        e.wall_p50_s *= 1.0 + f64::from(k) * 0.01;
        calib.record(&e, e.wall_p50_s * 1.1);
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
    gov.open_session(token(5, 1e9, 1e9)).expect("valid token");
    let gate = CancelGate::new();
    // Level 1: only the spill step fires.
    let l1 = gov
        .apply_memory_pressure(SessionId(5), 1, Some(&gate))
        .expect("session");
    assert_eq!(l1.len(), 1);
    assert_eq!(l1[0].step, DegradationStep::SpillColdArenas);
    assert!(!gate.is_requested(), "level 1 must not pause");
    // Level 3: all three fire IN THE DECLARED ORDER; pause requests the gate.
    let l3 = gov
        .apply_memory_pressure(SessionId(5), 3, Some(&gate))
        .expect("session");
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
    let events = gov.events();
    assert_eq!(events.len(), 4);
    assert!(events.windows(2).all(|w| w[0].ordinal < w[1].ordinal));
    assert!(events.iter().all(|e| !e.attribution.is_empty()));
    verdict(
        "ss-005",
        "ladder fires spill->coarsen->pause in declared order; snapshot round-trip exact",
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
fn ss_008_panicking_submission_releases_every_idempotency_waiter() {
    const WAITERS: usize = 8;
    let gov = Arc::new(Governor::new());
    gov.open_session(token(80, 1e9, 1e9)).expect("valid token");
    let executions = Arc::new(AtomicU32::new(0));
    let rendezvous = Arc::new(std::sync::Barrier::new(WAITERS + 1));
    let (started_tx, started_rx) = std::sync::mpsc::channel();
    let (done_tx, done_rx) = std::sync::mpsc::channel();
    let key = "panic-terminal".to_string();

    let owner = {
        let gov = Arc::clone(&gov);
        let executions = Arc::clone(&executions);
        let rendezvous = Arc::clone(&rendezvous);
        let done_tx = done_tx.clone();
        let key = key.clone();
        std::thread::spawn(move || {
            let outcome = gov.submit_once(SessionId(80), &key, || {
                executions.fetch_add(1, Ordering::SeqCst);
                started_tx.send(()).expect("test receiver alive");
                rendezvous.wait();
                panic!("seeded submission panic");
            });
            done_tx.send(outcome).expect("test receiver alive");
        })
    };

    started_rx
        .recv_timeout(std::time::Duration::from_secs(5))
        .expect("owner reached caller work after installing Pending");
    let mut waiters = Vec::new();
    for _ in 0..WAITERS {
        let gov = Arc::clone(&gov);
        let rendezvous = Arc::clone(&rendezvous);
        let done_tx = done_tx.clone();
        let key = key.clone();
        waiters.push(std::thread::spawn(move || {
            rendezvous.wait();
            let outcome = gov.submit_once(SessionId(80), &key, || {
                panic!("a duplicate must never execute");
            });
            done_tx.send(outcome).expect("test receiver alive");
        }));
    }
    drop(done_tx);

    let deadline = std::time::Instant::now() + std::time::Duration::from_secs(5);
    let mut failures = Vec::with_capacity(WAITERS + 1);
    for _ in 0..=WAITERS {
        let remaining = deadline.saturating_duration_since(std::time::Instant::now());
        let outcome = done_rx
            .recv_timeout(remaining)
            .expect("a panicking owner must not strand idempotency waiters")
            .expect("session remains valid");
        match outcome {
            SubmitOutcome::Failed { receipt, what } => failures.push((receipt, what)),
            other => panic!("every caller must observe the terminal failure, got {other:?}"),
        }
    }
    owner.join().expect("panic was contained by submit_once");
    for waiter in waiters {
        waiter.join().expect("waiter returned");
    }

    assert_eq!(executions.load(Ordering::SeqCst), 1, "work ran once");
    let (receipt, diagnosis) = &failures[0];
    assert!(
        failures.iter().all(|(r, d)| r == receipt && d == diagnosis),
        "owner and duplicates must share one failure receipt"
    );
    assert!(diagnosis.contains("seeded submission panic"), "{diagnosis}");
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
        matches!(retry, SubmitOutcome::Failed { receipt: r, .. } if r == *receipt),
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
        "seeded panic: one execution, one terminal failure receipt, all waiters released within 5s, no charge",
    );
}

#[test]
fn ss_009_invalid_resources_fail_closed_without_poisoning_meters() {
    let gov = Governor::new();
    let invalid = [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1.0];
    let mut next_id = 90u64;
    for value in invalid {
        for field in ["core_s", "wall_s", "cores"] {
            let id = next_id;
            next_id += 1;
            let mut bad = token(id, 100.0, 100.0);
            match field {
                "core_s" => bad.core_s = value,
                "wall_s" => bad.wall_s = value,
                "cores" => bad.cores = value,
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

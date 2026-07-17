//! Falsifier-pairing conformance (addendum Proposal 6, the qmao.4 bead).
//! Acceptance: no certificate class registers without ≥1 falsifier; the
//! consequence×doubt allocator concentrates budget correctly (monotone in
//! both factors, honest cold-start/floor boundaries); a reported discrepancy
//! produces escaped pending-candidate payloads; the catalog lint names
//! undeclared classes without pretending to authorize release; preliminary
//! class-level yield review is window-idempotent.

use fs_evidence::{
    AttemptRecord, ClaimContext, FalsifierAttempt, FalsifierHistory, FalsifierOutcome,
    FalsifierRegistry, FalsifierSpec, FalsifyError, allocate_budget,
    falsify::{
        CONSEQUENCE_FLOOR, DOUBT_COLD_START, DOUBT_FLOOR, MAX_FALSIFIERS_PER_CLASS, RENT_DECAY,
        RENT_SHARE_FLOOR, RENT_VOLUME,
    },
};

const SUITE: &str = "fs-evidence/falsify";
const FIXED_INPUT_SEED: u64 = 0;

fn verdict(case: &str, detail: &str) {
    let mut emitter = fs_obs::Emitter::new(SUITE, case);
    let event = emitter.emit(
        fs_obs::Severity::Info,
        fs_obs::EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: case.to_string(),
            pass: true,
            detail: detail.to_string(),
            seed: FIXED_INPUT_SEED,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("falsifier verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("falsifier verdict must use the fs-obs wire schema");
    println!("{line}");
}

fn claim(class: &str, regime: &str, consequence: f64) -> ClaimContext {
    ClaimContext {
        class: class.to_string(),
        regime: regime.to_string(),
        consequence,
    }
}

fn registry_for(entries: &[(&str, &str)]) -> FalsifierRegistry {
    let mut registry = FalsifierRegistry::new();
    for (class, falsifier) in entries {
        registry
            .register(
                class,
                vec![FalsifierSpec {
                    name: (*falsifier).to_string(),
                    method: "independent retained test checker".to_string(),
                }],
            )
            .expect("unique test declaration");
    }
    registry
}

fn attempt(
    attempt_id: impl Into<String>,
    class: &str,
    regime: &str,
    falsifier: &str,
    compute_s: f64,
    outcome: FalsifierOutcome,
) -> FalsifierAttempt {
    FalsifierAttempt {
        attempt_id: attempt_id.into(),
        class: class.to_string(),
        regime: regime.to_string(),
        falsifier: falsifier.to_string(),
        claim_revision: format!("claim-{class}-r1"),
        artifact_id: format!("artifact-{class}-r1"),
        seed: 7,
        compute_s,
        outcome,
    }
}

fn record_pass(
    history: &mut FalsifierHistory,
    registry: &FalsifierRegistry,
    attempt_id: impl Into<String>,
    class: &str,
    regime: &str,
    falsifier: &str,
    compute_s: f64,
) -> Result<AttemptRecord, FalsifyError> {
    history.record_attempt(
        registry,
        &attempt(
            attempt_id,
            class,
            regime,
            falsifier,
            compute_s,
            FalsifierOutcome::NoDiscrepancy,
        ),
    )
}

#[test]
fn fp_001_registration_requires_a_falsifier() {
    let mut r = FalsifierRegistry::new();
    // Empty falsifier list REFUSES — the rule at its source.
    let refusal = r.register("sampled-interface-agreement", Vec::new());
    assert_eq!(
        refusal,
        Err(FalsifyError::NoFalsifier {
            class: "sampled-interface-agreement".to_string()
        })
    );
    // With a falsifier it registers; duplicates refuse.
    r.register(
        "sampled-interface-agreement",
        vec![FalsifierSpec {
            name: "independent-retained-sample-replay".to_string(),
            method: "separately structured replay checker".to_string(),
        }],
    )
    .expect("registers");
    assert!(matches!(
        r.register(
            "sampled-interface-agreement",
            vec![FalsifierSpec {
                name: "x".to_string(),
                method: "y".to_string()
            }]
        ),
        Err(FalsifyError::Duplicate { .. })
    ));
    assert!(matches!(
        FalsifierRegistry::new().register(
            "x",
            vec![FalsifierSpec {
                name: " ".to_string(),
                method: "placeholder".to_string(),
            }]
        ),
        Err(FalsifyError::InvalidSpec { .. })
    ));
    assert!(matches!(
        FalsifierRegistry::new().register(
            "x",
            vec![
                FalsifierSpec {
                    name: "same".to_string(),
                    method: "method-a".to_string(),
                },
                FalsifierSpec {
                    name: "same".to_string(),
                    method: "method-b".to_string(),
                },
            ]
        ),
        Err(FalsifyError::DuplicateSpec { .. })
    ));
    let excessive = vec![
        FalsifierSpec {
            name: "bounded".to_string(),
            method: "bounded method".to_string(),
        };
        MAX_FALSIFIERS_PER_CLASS + 1
    ];
    assert!(matches!(
        FalsifierRegistry::new().register("too-many", excessive),
        Err(FalsifyError::ResourceLimit { .. })
    ));
    // The starting CATALOG covers both sampled and ambitious continuum claims.
    let std_reg = FalsifierRegistry::standard();
    for class in [
        "sampled-interface-agreement",
        "continuum-watertightness",
        "conservation",
        "adjoint-gradient",
        "surrogate-accept",
        "symmetry-block-solve",
        "validated-color",
    ] {
        let fs = std_reg.falsifiers(class).expect("registered");
        assert!(!fs.is_empty(), "{class} must carry a falsifier");
        assert!(!fs[0].method.is_empty(), "{class}: the method is stated");
    }
    assert!(std_reg.falsifiers("warp-drive").is_err());
    verdict(
        "fp-001",
        "empty/malformed/duplicate declarations refused; catalog separates sampled and continuum geometry claims",
    );
}

#[test]
#[allow(clippy::too_many_lines)] // One allocator conformance matrix shares one seeded history.
fn fp_002_budget_allocator_monotone_with_honest_boundaries() {
    let mut history = FalsifierHistory::new();
    let registry = registry_for(&[
        ("A", "test-check"),
        ("B", "test-check"),
        ("one", "test-check"),
    ]);
    // Build asymmetric history: class A has a strong record in regime r1,
    // class B has embarrassed itself.
    for index in 0..999 {
        record_pass(
            &mut history,
            &registry,
            format!("a-pass-{index}"),
            "A",
            "r1",
            "test-check",
            1.0,
        )
        .expect("pass");
    }
    for index in 0..79 {
        record_pass(
            &mut history,
            &registry,
            format!("b-pass-{index}"),
            "B",
            "r1",
            "test-check",
            1.0,
        )
        .expect("pass");
    }
    for hit_index in 0..20 {
        history
            .record_attempt(
                &registry,
                &attempt(
                    format!("b-hit-{hit_index}"),
                    "B",
                    "r1",
                    "test-check",
                    1.0,
                    FalsifierOutcome::Discrepancy {
                        detail: format!("flux mismatch {hit_index}"),
                    },
                ),
            )
            .expect("hit");
    }
    record_pass(
        &mut history,
        &registry,
        "one-pass-0",
        "one",
        "r1",
        "test-check",
        1.0,
    )
    .expect("pass");
    // Doubt: cold start is MAX; a large pass history reduces doubt gradually,
    // while one pass cannot collapse it to the floor.
    assert!(
        (history.doubt("C", "r1") - DOUBT_COLD_START).abs() < 1e-12,
        "cold start"
    );
    let doubt_a = history.doubt("A", "r1");
    assert!(
        doubt_a > DOUBT_FLOOR && doubt_a < 0.2,
        "large clean history remains conservatively above the floor: {doubt_a}"
    );
    assert!(
        history.doubt("one", "r1") > 0.9,
        "one pass must not manufacture confidence"
    );
    let doubt_b = history.doubt("B", "r1");
    assert!(doubt_b > 0.4, "an embarrassed class is doubted: {doubt_b}");
    // Monotone in doubt (same consequence).
    let claims = vec![
        claim("A", "r1", 5.0),
        claim("B", "r1", 5.0),
        claim("C", "r1", 5.0),
    ];
    let alloc = allocate_budget(100.0, &claims, &history).expect("allocation");
    assert!(
        (alloc.iter().sum::<f64>() - 100.0).abs() < 1e-9,
        "budget fully spent"
    );
    assert!(
        alloc[0] < alloc[1] && alloc[1] < alloc[2],
        "doubt-monotone: {alloc:?} (A trusted < B embarrassed < C unknown)"
    );
    // Monotone in consequence (same class/doubt).
    let claims2 = vec![claim("B", "r1", 1.0), claim("B", "r1", 10.0)];
    let alloc2 = allocate_budget(50.0, &claims2, &history).expect("allocation");
    assert!(alloc2[0] < alloc2[1], "consequence-monotone: {alloc2:?}");
    assert!(
        (alloc2[1] / alloc2[0] - 10.0).abs() < 1e-9,
        "proportional to consequence"
    );
    // Boundaries: zero-dependents claim gets a nonzero floor share;
    // zero claims spend zero.
    let claims3 = vec![claim("B", "r1", 0.0), claim("B", "r1", 1.0)];
    let alloc3 = allocate_budget(10.0, &claims3, &history).expect("allocation");
    assert!(alloc3[0] > 0.0, "no-dependents claim keeps a floor share");
    assert!(
        (alloc3[0] / alloc3[1] - CONSEQUENCE_FLOOR).abs() < 1e-9,
        "floor ratio is the declared constant"
    );
    assert!(
        allocate_budget(10.0, &[], &history)
            .expect("empty allocation")
            .is_empty(),
        "no claims, no spend"
    );
    assert!(allocate_budget(f64::NAN, &claims, &history).is_err());
    assert!(allocate_budget(f64::INFINITY, &claims, &history).is_err());
    assert!(allocate_budget(1.0, &[claim("A", "r1", f64::MAX)], &history).is_ok());
    let extremes = allocate_budget(
        1.0,
        &[
            claim("extreme-a", "r1", f64::MAX),
            claim("extreme-b", "r1", f64::MAX),
        ],
        &FalsifierHistory::new(),
    )
    .expect("max-rescaled extreme allocation");
    assert!(extremes.iter().all(|value| value.is_finite()));
    assert!((extremes.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    assert!(
        allocate_budget(0.0, &[claim("bad id", "r1", f64::NAN)], &history).is_err(),
        "zero budget must not bypass claim validation"
    );
    verdict(
        "fp-002",
        "monotone in consequence and doubt; cold-start max, asymptotic doubt floor, \
         dependent-free floor, empty-job zero",
    );
}

#[test]
#[allow(clippy::too_many_lines)] // One discrepancy lifecycle shares one idempotent history.
fn fp_003_hit_wiring_creates_tombstone_and_bug_report() {
    let mut history = FalsifierHistory::new();
    let registry = FalsifierRegistry::standard();
    let hit = attempt(
        "attempt-adjoint-7",
        "adjoint-gradient",
        "Re~1e3",
        "finite-difference-spot-check",
        2.5,
        FalsifierOutcome::Discrepancy {
            detail: "FD says 0.031, tape says \"0.29\"\nnext".to_string(),
        },
    );
    let record = history.record_attempt(&registry, &hit).expect("valid hit");
    assert!(!record.duplicate);
    let tombstone = record.tombstone.expect("candidate tombstone");
    let bug = record.estimator_bug.expect("candidate bug");
    for (name, json) in [("tombstone", tombstone.json()), ("bug", bug.json())] {
        assert!(json.contains("0.031"), "{name} carries the evidence");
        assert!(!json.contains("\nnext"), "{name} contains no raw newline");
        assert!(json.contains("\\nnext"), "{name} escapes the newline");
        assert!(json.contains("\"schema_id\":\"fs-evidence/falsifier-candidate\""));
        assert!(json.contains("\"schema_version\":1"));
    }
    assert_eq!(
        tombstone.json(),
        "{\"schema_id\":\"fs-evidence/falsifier-candidate\",\"schema_version\":1,\"kind\":\"tombstone-candidate\",\"source\":\"falsifier-attempt-candidate\",\"adjudication\":\"pending\",\"attempt_id\":\"attempt-adjoint-7\",\"class\":\"adjoint-gradient\",\"regime\":\"Re~1e3\",\"falsifier\":\"finite-difference-spot-check\",\"claim_revision\":\"claim-adjoint-gradient-r1\",\"artifact_id\":\"artifact-adjoint-gradient-r1\",\"seed_bits\":\"0000000000000007\",\"compute_s_bits\":\"4004000000000000\",\"detail\":\"FD says 0.031, tape says \\\"0.29\\\"\\nnext\"}"
    );
    assert_eq!(
        bug.json(),
        "{\"schema_id\":\"fs-evidence/falsifier-candidate\",\"schema_version\":1,\"kind\":\"estimator-bug-candidate\",\"adjudication\":\"pending\",\"attempt_id\":\"attempt-adjoint-7\",\"class\":\"adjoint-gradient\",\"regime\":\"Re~1e3\",\"caught_by\":\"finite-difference-spot-check\",\"claim_revision\":\"claim-adjoint-gradient-r1\",\"artifact_id\":\"artifact-adjoint-gradient-r1\",\"seed_bits\":\"0000000000000007\",\"compute_s_bits\":\"4004000000000000\",\"evidence\":\"FD says 0.031, tape says \\\"0.29\\\"\\nnext\"}"
    );
    assert!(
        tombstone
            .json()
            .contains("\"kind\":\"tombstone-candidate\"")
    );
    assert!(bug.json().contains("\"kind\":\"estimator-bug-candidate\""));
    assert!(tombstone.json().contains("\"adjudication\":\"pending\""));
    assert!(bug.json().contains("\"adjudication\":\"pending\""));
    let replay = history
        .record_attempt(&registry, &hit)
        .expect("identical retry is idempotent");
    assert!(replay.duplicate);
    assert_eq!(replay.tombstone.expect("replayed candidate"), tombstone);
    // The hit moved the class's doubt and yield.
    assert!(history.doubt(&hit.class, &hit.regime) > 0.9);
    let (hits, compute, runs) = history.yield_of(&hit.class).expect("yield");
    assert_eq!((hits, runs), (1, 1));
    assert!((compute - 2.5).abs() < 1e-12);
    let mut invalid_compute = hit.clone();
    invalid_compute.attempt_id = "attempt-negative".to_string();
    invalid_compute.compute_s = -1.0;
    assert!(history.record_attempt(&registry, &invalid_compute).is_err());
    invalid_compute.attempt_id = "attempt-nan".to_string();
    invalid_compute.compute_s = f64::NAN;
    assert!(history.record_attempt(&registry, &invalid_compute).is_err());
    invalid_compute.attempt_id = "attempt-zero".to_string();
    invalid_compute.compute_s = 0.0;
    assert!(history.record_attempt(&registry, &invalid_compute).is_err());
    let mut invalid_id = hit.clone();
    invalid_id.attempt_id = "attempt-bad-id".to_string();
    invalid_id.class = "bad id".to_string();
    assert!(history.record_attempt(&registry, &invalid_id).is_err());
    let mut empty_detail = hit.clone();
    empty_detail.attempt_id = "attempt-empty-detail".to_string();
    empty_detail.outcome = FalsifierOutcome::Discrepancy {
        detail: "  ".to_string(),
    };
    assert!(matches!(
        history.record_attempt(&registry, &empty_detail),
        Err(FalsifyError::InvalidText { .. })
    ));
    let mut collision = hit.clone();
    collision.compute_s = 3.0;
    assert!(matches!(
        history.record_attempt(&registry, &collision),
        Err(FalsifyError::AttemptCollision { .. })
    ));
    let (hits_after_replay, compute_after_replay, runs_after_replay) =
        history.yield_of(&hit.class).expect("yield after retry");
    assert_eq!((hits_after_replay, runs_after_replay), (1, 1));
    assert!((compute_after_replay - 2.5).abs() < 1e-12);

    let replaced_catalog = FalsifierRegistry::new();
    let replay_after_catalog_replacement = history
        .record_attempt(&replaced_catalog, &hit)
        .expect("an accepted exact retry is independent of the current catalog");
    assert!(replay_after_catalog_replacement.duplicate);
    assert_eq!(
        history
            .yield_of(&hit.class)
            .expect("yield after catalog replacement"),
        (1, 2.5, 1),
        "an idempotent retry must not mutate diagnostic yield"
    );
    let mut replaced_catalog_collision = hit.clone();
    replaced_catalog_collision.compute_s = 4.0;
    assert!(matches!(
        history.record_attempt(&replaced_catalog, &replaced_catalog_collision),
        Err(FalsifyError::AttemptCollision { .. })
    ));
    verdict(
        "fp-003",
        "reported discrepancy -> escaped pending candidates; diagnostic doubt/yield updated",
    );
}

#[test]
fn fp_004_catalog_lint_is_not_release_authority() {
    let registry = FalsifierRegistry::standard();
    // Declared classes satisfy catalog completeness only.
    assert!(
        registry
            .catalog_gate(&["sampled-interface-agreement", "conservation"])
            .expect("bounded valid catalog query")
            .is_empty(),
        "declared classes pass the catalog lint"
    );
    let violations = registry
        .catalog_gate(&["sampled-interface-agreement", "novel-certificate-v0"])
        .expect("bounded valid catalog query");
    assert_eq!(violations, vec!["novel-certificate-v0".to_string()]);
    assert_eq!(
        registry
            .catalog_gate(&["z-missing", "a-missing", "z-missing"])
            .expect("duplicate bounded catalog query"),
        vec!["a-missing".to_string(), "z-missing".to_string()],
        "catalog diagnostics are canonical unique classes, not an amplified caller multiset"
    );
    assert!(matches!(
        registry.catalog_gate(&[" "]),
        Err(FalsifyError::InvalidIdentifier { .. })
    ));
    verdict(
        "fp-004",
        "catalog lint names undeclared classes but makes no exact-instance or release claim",
    );
}

#[test]
#[allow(clippy::too_many_lines)] // One stateful window sequence proves decay and idempotence together.
fn fp_005_rent_review_decays_yieldless_falsifiers() {
    let mut history = FalsifierHistory::new();
    let registry = registry_for(&[
        ("quiet", "test-check"),
        ("worker", "test-check"),
        ("young", "test-check"),
    ]);
    // "quiet" runs at meaningful volume with zero hits; "worker" reports a
    // candidate discrepancy in its first governance window.
    for index in 0..RENT_VOLUME {
        record_pass(
            &mut history,
            &registry,
            format!("quiet-initial-{index}"),
            "quiet",
            "r",
            "test-check",
            0.5,
        )
        .expect("pass");
    }
    for index in 0..RENT_VOLUME - 1 {
        record_pass(
            &mut history,
            &registry,
            format!("worker-initial-{index}"),
            "worker",
            "r",
            "test-check",
            0.5,
        )
        .expect("pass");
    }
    history
        .record_attempt(
            &registry,
            &attempt(
                "worker-initial-hit",
                "worker",
                "r",
                "test-check",
                0.5,
                FalsifierOutcome::Discrepancy {
                    detail: "reported one discrepancy".to_string(),
                },
            ),
        )
        .expect("hit");
    // Below-volume class is exempt.
    for index in 0..10 {
        record_pass(
            &mut history,
            &registry,
            format!("young-{index}"),
            "young",
            "r",
            "test-check",
            0.5,
        )
        .expect("pass");
    }
    let decayed = history.rent_review().expect("review");
    assert_eq!(
        decayed.len(),
        1,
        "only the yield-less at-volume class decays: {decayed:?}"
    );
    assert_eq!(decayed[0].0, "quiet");
    assert!((history.share("quiet") - 0.5).abs() < 1e-12);
    assert!(
        (history.share("worker") - 1.0).abs() < 1e-12,
        "classes with a reported discrepancy keep their first-window share"
    );
    assert!(
        (history.share("young") - 1.0).abs() < 1e-12,
        "low-volume classes are exempt"
    );
    assert!(
        history.rent_review().expect("idempotent review").is_empty(),
        "unchanged repeated review must not decay again"
    );
    let claims = vec![claim("quiet", "r", 1.0), claim("worker", "r", 1.0)];
    let alloc = allocate_budget(10.0, &claims, &history).expect("allocation");
    assert!(
        alloc[0] < alloc[1],
        "decayed share reduces allocation: {alloc:?}"
    );
    // Every new zero-hit window may decay again to the floor, while the old
    // worker hit does not grant permanent immunity.
    for window in 0..10 {
        for index in 0..RENT_VOLUME {
            record_pass(
                &mut history,
                &registry,
                format!("quiet-window-{window}-{index}"),
                "quiet",
                "r",
                "test-check",
                0.5,
            )
            .expect("pass");
            record_pass(
                &mut history,
                &registry,
                format!("worker-window-{window}-{index}"),
                "worker",
                "r",
                "test-check",
                0.5,
            )
            .expect("pass");
        }
        let _ = history.rent_review().expect("new-window review");
    }
    assert!(
        history.share("quiet") >= RENT_SHARE_FLOOR,
        "share floors at {RENT_SHARE_FLOOR}, never zero"
    );
    assert!(
        history.share("worker") < 1.0,
        "one historical hit does not immunize future zero-hit windows"
    );
    verdict(
        "fp-005",
        "class-level review is window-idempotent; new zero-hit windows decay to a floor; old hits do not immunize",
    );
}

#[test]
fn fp_006_subthreshold_review_cannot_erase_the_window() {
    let mut history = FalsifierHistory::new();
    let registry = registry_for(&[("drip", "test-check")]);
    for index in 0..RENT_VOLUME - 1 {
        record_pass(
            &mut history,
            &registry,
            format!("drip-{index}"),
            "drip",
            "r",
            "test-check",
            0.5,
        )
        .expect("pass");
    }
    assert!(
        history
            .rent_review()
            .expect("subthreshold review")
            .is_empty()
    );
    assert_eq!(history.share("drip").to_bits(), 1.0_f64.to_bits());
    record_pass(
        &mut history,
        &registry,
        "drip-99",
        "drip",
        "r",
        "test-check",
        0.5,
    )
    .expect("threshold pass");
    let decayed = history.rent_review().expect("threshold review");
    assert_eq!(decayed, vec![("drip".to_string(), 0.5)]);
    verdict(
        "fp-006",
        "reviewing at 99 runs cannot erase them; the 100th closes the diagnostic window",
    );
}

#[test]
fn fp_007_rent_windows_are_review_schedule_invariant() {
    let registry = registry_for(&[("batch", "test-check")]);
    let mut batched = FalsifierHistory::new();
    let mut incremental = FalsifierHistory::new();
    for index in 0..(2 * RENT_VOLUME) {
        record_pass(
            &mut batched,
            &registry,
            format!("batched-{index}"),
            "batch",
            "r",
            "test-check",
            0.5,
        )
        .expect("batched pass");
        record_pass(
            &mut incremental,
            &registry,
            format!("incremental-{index}"),
            "batch",
            "r",
            "test-check",
            0.5,
        )
        .expect("incremental pass");
        if index + 1 == RENT_VOLUME {
            let _ = incremental.rent_review().expect("first fixed window");
        }
    }
    let batched_decays = batched.rent_review().expect("two batched windows");
    let incremental_decays = incremental.rent_review().expect("second fixed window");
    assert_eq!(batched.share("batch").to_bits(), 0.25_f64.to_bits());
    assert_eq!(incremental.share("batch").to_bits(), 0.25_f64.to_bits());
    assert_eq!(batched_decays, vec![("batch".to_string(), 0.25)]);
    assert_eq!(incremental_decays, vec![("batch".to_string(), 0.25)]);
    verdict(
        "fp-007",
        "two fixed-volume clean windows decay identically under batched and incremental review schedules",
    );
}

#[test]
#[allow(clippy::too_many_lines)] // One atomicity matrix must compare every mutation to one baseline.
fn fp_008_attempt_identity_is_atomic_and_full_width() {
    let mut registry = FalsifierRegistry::new();
    registry
        .register(
            "class-a",
            vec![
                FalsifierSpec {
                    name: "check-a".to_string(),
                    method: "independent retained checker A".to_string(),
                },
                FalsifierSpec {
                    name: "check-alt".to_string(),
                    method: "independent retained checker alternate".to_string(),
                },
            ],
        )
        .expect("class A");
    registry
        .register(
            "class-b",
            vec![FalsifierSpec {
                name: "check-a".to_string(),
                method: "independent retained checker B".to_string(),
            }],
        )
        .expect("class B");

    let mut history = FalsifierHistory::new();
    let baseline = attempt(
        "immutable-attempt",
        "class-a",
        "regime-a",
        "check-a",
        1.0,
        FalsifierOutcome::NoDiscrepancy,
    );
    history
        .record_attempt(&registry, &baseline)
        .expect("baseline attempt");
    let stable_state = format!("{history:?}");

    let mut mutations = Vec::new();
    let mut changed = baseline.clone();
    changed.class = "class-b".to_string();
    mutations.push(changed);
    let mut changed = baseline.clone();
    changed.regime = "regime-b".to_string();
    mutations.push(changed);
    let mut changed = baseline.clone();
    changed.falsifier = "check-alt".to_string();
    mutations.push(changed);
    let mut changed = baseline.clone();
    changed.claim_revision = "claim-class-a-r2".to_string();
    mutations.push(changed);
    let mut changed = baseline.clone();
    changed.artifact_id = "artifact-class-a-r2".to_string();
    mutations.push(changed);
    let mut changed = baseline.clone();
    changed.seed = u64::MAX;
    mutations.push(changed);
    let mut changed = baseline.clone();
    changed.compute_s = 2.0;
    mutations.push(changed);
    let mut changed = baseline.clone();
    changed.outcome = FalsifierOutcome::Discrepancy {
        detail: "different outcome".to_string(),
    };
    mutations.push(changed);

    for changed in mutations {
        assert!(matches!(
            history.record_attempt(&registry, &changed),
            Err(FalsifyError::AttemptCollision { .. })
        ));
        assert_eq!(format!("{history:?}"), stable_state);
    }

    let mut unknown = baseline.clone();
    unknown.attempt_id = "unknown-class-attempt".to_string();
    unknown.class = "class-unknown".to_string();
    assert!(matches!(
        history.record_attempt(&registry, &unknown),
        Err(FalsifyError::Unknown { .. })
    ));
    assert_eq!(format!("{history:?}"), stable_state);

    let mut undeclared = baseline.clone();
    undeclared.attempt_id = "undeclared-falsifier-attempt".to_string();
    undeclared.falsifier = "check-missing".to_string();
    assert!(matches!(
        history.record_attempt(&registry, &undeclared),
        Err(FalsifyError::UnregisteredFalsifier { .. })
    ));
    assert_eq!(format!("{history:?}"), stable_state);

    let mut full_width = attempt(
        "full-width-seed",
        "class-a",
        "regime-a",
        "check-a",
        1.0,
        FalsifierOutcome::Discrepancy {
            detail: "full-width seed transport".to_string(),
        },
    );
    full_width.seed = u64::MAX;
    let candidate = history
        .record_attempt(&registry, &full_width)
        .expect("full-width attempt")
        .tombstone
        .expect("candidate");
    assert!(
        candidate
            .json()
            .contains("\"seed_bits\":\"ffffffffffffffff\"")
    );

    let mut boundary_history = FalsifierHistory::new();
    for index in 0..RENT_VOLUME - 1 {
        record_pass(
            &mut boundary_history,
            &registry,
            format!("boundary-{index}"),
            "class-a",
            "regime-a",
            "check-a",
            1.0,
        )
        .expect("pre-boundary pass");
    }
    let boundary = attempt(
        "boundary-final",
        "class-a",
        "regime-a",
        "check-a",
        1.0,
        FalsifierOutcome::NoDiscrepancy,
    );
    assert!(
        !boundary_history
            .record_attempt(&registry, &boundary)
            .expect("boundary close")
            .duplicate
    );
    assert!(
        boundary_history
            .record_attempt(&registry, &boundary)
            .expect("idempotent boundary retry")
            .duplicate
    );
    assert_eq!(
        boundary_history.rent_review().expect("one closed window"),
        vec![("class-a".to_string(), RENT_DECAY)]
    );
    assert!(
        boundary_history
            .rent_review()
            .expect("no duplicate window")
            .is_empty()
    );

    verdict(
        "fp-008",
        "all attempt-identity fields are collision-atomic; unknown/undeclared inputs roll back; full-width seed and boundary replay are lossless",
    );
}

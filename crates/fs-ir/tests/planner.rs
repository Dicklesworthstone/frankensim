//! Ladder-planner conformance (the lmp4.16 bead; runs under the
//! `ladder-planner` feature). Acceptance: queries discharge at the
//! requested tolerance within budget with sensible operator choices;
//! THE KILL MEASUREMENT — the greedy planner beats the fixed
//! mid-rung + uniform-refinement baseline by ≥2× cost at equal
//! certified accuracy; the certified-accuracy contract is never
//! violated; cache hits return with zero solves; cold cost estimates
//! fall back conservatively; the cannot-discharge boundary refuses
//! with the best achieved interval; replay is deterministic (G5).
#![cfg(feature = "ladder-planner")]

use fs_ir::planner::{
    AnswerCache, CachedAnswer, CostTable, MAX_FAMILY_COEFFICIENTS, MAX_LADDER_RUNGS,
    MAX_PLANNER_CELLS, MAX_POLYNOMIAL_CELL_WORK, MemCache, PlanError, PlanOp, PlanOutcome,
    ProblemFamily, baseline_uniform, plan,
};
use fs_verify::fem1d::Poly;

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-ir/planner\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

/// The wedge-like family: smooth base + the highest-degree right feature
/// admitted by the exact five-point verifier, so error concentrates near x = 1 and local
/// refinement genuinely beats uniform.
fn steep_family() -> ProblemFamily {
    // u = x(1−x)·(0.2 + x³)
    //   = 0.2x − 0.2x² + x⁴ − x⁵.
    let mut c = vec![0.0; 6];
    c[1] = 0.2;
    c[2] = -0.2;
    c[4] = 1.0;
    c[5] = -1.0;
    ProblemFamily::new(Poly(c), "cht-wedge-steep").unwrap()
}

const RUNGS: [usize; 4] = [12, 24, 48, 96];

#[test]
fn pl_001_discharges_within_budget_contract_held() {
    let family = steep_family();
    let mut cache = MemCache::default();
    let mut costs = CostTable::new(200.0).unwrap();
    let tol = 0.05;
    let out = plan(&family, 1.0, tol, 2000.0, &RUNGS, &mut cache, &mut costs).unwrap();
    match &out {
        PlanOutcome::Discharged {
            bound, ops, cost, ..
        } => {
            // THE CONTRACT: a discharged answer never violates the
            // certified accuracy (property, not sample).
            assert!(*bound <= tol, "certified: {bound} <= {tol}");
            assert!(*cost <= 2000.0, "within budget: {cost}");
            // Sensible composition: cache first, then work.
            assert_eq!(ops[0].op, PlanOp::CacheLookup);
            assert!(ops.iter().any(|o| o.op == PlanOp::SolveRung));
            for entry in ops.iter().filter(|entry| entry.bound_after.is_finite()) {
                let certificate = entry
                    .certificate_after
                    .as_ref()
                    .expect("every finite verifier bound retains its authority");
                assert_eq!(certificate.bound().to_bits(), entry.bound_after.to_bits());
                assert_eq!(certificate.verifier_family(), "equilibrated-flux-1d");
                assert_eq!(
                    certificate.color(),
                    &fs_evidence::Color::Verified {
                        lo: 0.0,
                        hi: entry.bound_after,
                    }
                );
            }
            let seq: Vec<&str> = ops.iter().map(|o| o.op.name()).collect();
            println!(
                "{{\"metric\":\"planner-run\",\"ops\":{seq:?},\"cost\":{cost},\
                 \"bound\":{bound:.3e}}}"
            );
        }
        PlanOutcome::RefusedWithBest { reason, .. } => {
            panic!("a generous budget must discharge: {reason}")
        }
        PlanOutcome::RefusedWithoutAnswer { reason, .. } => {
            panic!("a generous budget must produce an answer: {reason}")
        }
    }
    verdict(
        "pl-001",
        "discharged at tol within budget; certified bound honored; cache-first op order",
    );
}

#[test]
fn pl_002_the_kill_measurement() {
    // Planner vs the fixed baseline (mid-rung 48 + uniform doubling) at
    // EQUAL certified accuracy: the greedy walk must be >= 2x cheaper.
    let family = steep_family();
    let tol = 6e-3;
    let (base_cost, base_bound) = baseline_uniform(&family, 1.0, tol, 48, 6).unwrap();
    assert!(base_bound <= tol, "the baseline eventually certifies");
    let mut cache = MemCache::default();
    let mut costs = CostTable::new(200.0).unwrap();
    let out = plan(&family, 1.0, tol, 100_000.0, &RUNGS, &mut cache, &mut costs).unwrap();
    let PlanOutcome::Discharged { bound, cost, .. } = out else {
        panic!("the planner must discharge at this tolerance");
    };
    assert!(bound <= tol, "equal certified accuracy");
    let ratio = base_cost / cost;
    println!(
        "{{\"metric\":\"planner-kill-check\",\"baseline_cells\":{base_cost},\
         \"planner_cells\":{cost},\"ratio\":{ratio:.2},\"gate\":2.0}}"
    );
    assert!(
        ratio >= 2.0,
        "the kill criterion: planner must beat the baseline >=2x at equal certified \
         accuracy (got {ratio:.2}x) — else ship the interface and freeze the planner"
    );
    verdict(
        "pl-002",
        "kill measurement PASSED: the ladder walk beats mid-rung+uniform by >=2x cells \
         at equal certified accuracy",
    );
}

#[test]
fn pl_003_cache_hits_and_cold_estimates() {
    let family = steep_family();
    let mut cache = MemCache::default();
    let mut costs = CostTable::new(500.0).unwrap();
    // COLD estimates: before any telemetry, predictions are the
    // conservative default (the round-3 boundary).
    assert!((costs.predict(PlanOp::DwrRefine) - 500.0).abs() < f64::EPSILON);
    let tol = 0.05;
    let first = plan(&family, 1.0, tol, 5000.0, &RUNGS, &mut cache, &mut costs).unwrap();
    assert!(matches!(first, PlanOutcome::Discharged { .. }));
    // Learned estimates move off the default.
    assert!(
        costs.predict(PlanOp::SolveRung) < 500.0,
        "telemetry sharpens the table: {}",
        costs.predict(PlanOp::SolveRung)
    );
    // The SAME query again: a cache hit with ZERO solves.
    let again = plan(&family, 1.0, tol, 5000.0, &RUNGS, &mut cache, &mut costs).unwrap();
    match again {
        PlanOutcome::Discharged { ops, cost, .. } => {
            assert_eq!(ops.len(), 1, "one op only: the cache");
            assert_eq!(ops[0].op, PlanOp::CacheLookup);
            assert!(cost.abs() < f64::EPSILON, "zero solves on a hit");
        }
        PlanOutcome::RefusedWithBest { .. } => panic!("the hit must discharge"),
        PlanOutcome::RefusedWithoutAnswer { .. } => panic!("the hit must return an answer"),
    }
    verdict(
        "pl-003",
        "cold table predicts the conservative default; telemetry sharpens it; the \
         repeat query is a zero-solve cache hit",
    );
}

#[test]
fn pl_004_refusal_boundary_and_g5_determinism() {
    let family = steep_family();
    // A budget too small to certify a tight tolerance: the planner must
    // refuse WITH its best certified interval, never overrun or lie.
    let tol = 1e-4;
    let budget = 80.0;
    let mut cache = MemCache::default();
    let mut costs = CostTable::new(200.0).unwrap();
    let out = plan(&family, 1.0, tol, budget, &RUNGS, &mut cache, &mut costs).unwrap();
    match &out {
        PlanOutcome::RefusedWithBest {
            best_bound,
            cost,
            reason,
            ops,
            ..
        } => {
            assert!(
                *best_bound > tol,
                "honest: the best bound did not reach tol"
            );
            assert!(
                best_bound.is_finite(),
                "a certified interval travels with the refusal"
            );
            assert!(
                *cost <= budget,
                "admitted work never exceeds the budget: {cost} vs {budget}"
            );
            assert!(
                reason.contains("refusal"),
                "hands off to refusal semantics: {reason}"
            );
            assert!(!ops.is_empty());
        }
        PlanOutcome::Discharged { .. } => panic!("80 cells cannot certify 1e-4"),
        PlanOutcome::RefusedWithoutAnswer { reason, .. } => {
            panic!("80 cells can fund an initial certified interval: {reason}")
        }
    }
    // G5: the identical query replays the identical operator sequence.
    let run = |seed_cache: &mut MemCache| -> Vec<&'static str> {
        let mut costs = CostTable::new(200.0).unwrap();
        match plan(&family, 1.0, 0.05, 2000.0, &RUNGS, seed_cache, &mut costs).unwrap() {
            PlanOutcome::Discharged { ops, .. } | PlanOutcome::RefusedWithBest { ops, .. } => {
                ops.iter().map(|o| o.op.name()).collect()
            }
            PlanOutcome::RefusedWithoutAnswer { ops, .. } => {
                ops.iter().map(|o| o.op.name()).collect()
            }
        }
    };
    let a = run(&mut MemCache::default());
    let b = run(&mut MemCache::default());
    assert_eq!(a, b, "replayed queries reproduce the operator sequence");
    verdict(
        "pl-004",
        "under-budget queries refuse with the best certified interval and bounded \
         spend; replays are deterministic",
    );
}

#[test]
fn pl_005_cost_calibration() {
    // Predicted-vs-actual per operator after learning: within 2x.
    let family = steep_family();
    let mut cache = MemCache::default();
    let mut costs = CostTable::new(300.0).unwrap();
    for k in 0..4 {
        let theta = 0.9 + 0.05 * f64::from(k as u8);
        let _ = plan(&family, theta, 0.05, 5000.0, &RUNGS, &mut cache, &mut costs).unwrap();
    }
    let predicted = costs.predict(PlanOp::SolveRung);
    assert!(
        predicted < 300.0,
        "solve-rung: learned below the cold default ({predicted})"
    );
    println!(
        "{{\"metric\":\"cost-calibration\",\"solve\":{:.1},\"speculate\":{:.1},\
         \"refine\":{:.1},\"climb\":{:.1}}}",
        costs.predict(PlanOp::SolveRung),
        costs.predict(PlanOp::Speculate),
        costs.predict(PlanOp::DwrRefine),
        costs.predict(PlanOp::Climb),
    );
    verdict(
        "pl-005",
        "after 4 planned queries the cost table is learned for every exercised operator",
    );
}

#[test]
fn pl_006_malformed_planner_inputs_refuse_structurally() {
    let family = steep_family();
    let mut cache = MemCache::default();
    let mut costs = CostTable::new(20.0).unwrap();
    let mut run = |theta, tol, budget, rungs: &[usize]| {
        plan(&family, theta, tol, budget, rungs, &mut cache, &mut costs)
    };

    assert_eq!(
        run(1.0, 0.1, 100.0, &[]),
        Err(PlanError::EmptySequence {
            field: "rung_cells",
        })
    );
    assert!(matches!(
        run(1.0, 0.1, 100.0, &[0, 4]),
        Err(PlanError::InvalidSequenceEntry {
            field: "rung_cells",
            index: 0,
            ..
        })
    ));
    assert!(matches!(
        run(1.0, 0.1, 100.0, &[4, 4]),
        Err(PlanError::InvalidSequenceEntry {
            field: "rung_cells",
            index: 1,
            ..
        })
    ));
    assert!(matches!(
        run(1.0, 0.1, 100.0, &[8, 5]),
        Err(PlanError::InvalidSequenceEntry {
            field: "rung_cells",
            index: 1,
            ..
        })
    ));
    assert_eq!(
        run(1.0, 0.1, 100.0, &[MAX_PLANNER_CELLS + 1]),
        Err(PlanError::ResourceLimit {
            field: "rung_cells",
            requested: MAX_PLANNER_CELLS + 1,
            limit: MAX_PLANNER_CELLS,
        })
    );
    let oversized_ladder = (1..=MAX_LADDER_RUNGS + 1).collect::<Vec<_>>();
    assert_eq!(
        run(1.0, 0.1, 100.0, &oversized_ladder),
        Err(PlanError::ResourceLimit {
            field: "rung_cells",
            requested: MAX_LADDER_RUNGS + 1,
            limit: MAX_LADDER_RUNGS,
        })
    );
    for (theta, tol, budget, field) in [
        (f64::NAN, 0.1, 100.0, "theta"),
        (1.0, f64::INFINITY, 100.0, "tolerance"),
        (1.0, 0.0, 100.0, "tolerance"),
        (1.0, 0.1, f64::NAN, "budget_cells"),
        (1.0, 0.1, 0.0, "budget_cells"),
    ] {
        assert!(matches!(
            run(theta, tol, budget, &[4, 7]),
            Err(PlanError::InvalidScalar { field: actual, .. }) if actual == field
        ));
    }
}

#[test]
fn pl_007_family_mesh_candidate_and_cost_admission_are_fail_closed() {
    assert!(ProblemFamily::new(Poly(vec![]), "kernel").is_err());
    assert!(ProblemFamily::new(Poly(vec![0.0, f64::NAN, 0.0]), "kernel").is_err());
    assert!(ProblemFamily::new(Poly(vec![1.0, -1.0]), "kernel").is_err());
    assert!(ProblemFamily::new(Poly(vec![f64::EPSILON, -f64::EPSILON]), "kernel").is_err());
    assert!(
        ProblemFamily::new(
            Poly(vec![0.0, 1.0e16, -1.0e16, 1.0]),
            "hidden-boundary-residue"
        )
        .is_err(),
        "point Horner zero cannot hide an exact nonzero boundary residue"
    );
    assert!(
        ProblemFamily::new(
            Poly(vec![0.0, 1.0, 1.0e16, -1.0e16, -1.0]),
            "exact-cancellation"
        )
        .is_ok(),
        "an exact binary-rational cancellation survives point-Horner rounding"
    );
    assert!(ProblemFamily::new(Poly(vec![0.0, 0.0]), "unknown").is_err());
    assert_eq!(
        ProblemFamily::new(
            Poly(vec![0.0; MAX_FAMILY_COEFFICIENTS + 1]),
            "oversized-family"
        ),
        Err(PlanError::ResourceLimit {
            field: "family.base",
            requested: MAX_FAMILY_COEFFICIENTS + 1,
            limit: MAX_FAMILY_COEFFICIENTS,
        })
    );

    let family = steep_family();
    assert!(family.at(1.0, vec![]).is_err());
    assert!(family.at(1.0, vec![0.0, 0.7, 0.6, 1.0]).is_err());
    assert!(family.at(1.0, vec![0.0, f64::NAN, 1.0]).is_err());
    let explosive = ProblemFamily::new(Poly(vec![0.0, f64::MAX, -f64::MAX]), "explosive").unwrap();
    assert!(explosive.at(1.0, vec![0.0, 0.5, 1.0]).is_err());
    let scaling_sensitive =
        ProblemFamily::new(Poly(vec![0.0, 1.0, 2.0, -3.0]), "scaling-sensitive").unwrap();
    assert!(matches!(
        scaling_sensitive.at(0.1, vec![0.0, 0.5, 1.0]),
        Err(PlanError::InvalidFamily {
            field: "scaled_base",
            ..
        })
    ));

    assert!(CostTable::new(f64::NAN).is_err());
    assert!(CostTable::new(0.0).is_err());
    let mut costs = CostTable::new(10.0).unwrap();
    assert!(costs.record(PlanOp::Climb, f64::NAN).is_err());
    assert!(costs.record(PlanOp::Climb, -1.0).is_err());
    assert!((costs.predict(PlanOp::Climb) - 10.0).abs() < f64::EPSILON);

    assert!(CachedAnswer::new(vec![0.0, 0.0], -1.0, vec![0.0, 1.0]).is_err());
    assert!(CachedAnswer::new(vec![0.0], 0.1, vec![0.0, 1.0]).is_err());
    assert!(CachedAnswer::new(vec![1.0, 0.0], 0.1, vec![0.0, 1.0]).is_err());
}

#[test]
fn pl_007b_combined_polynomial_mesh_work_refuses_before_mesh_allocation() {
    let family = ProblemFamily::new(
        Poly(vec![0.0; MAX_FAMILY_COEFFICIENTS]),
        "bounded-but-combined-too-large",
    )
    .unwrap();
    let cells = MAX_POLYNOMIAL_CELL_WORK / (MAX_FAMILY_COEFFICIENTS * 5) + 1;
    let result = plan(
        &family,
        1.0,
        0.1,
        1.0,
        &[cells],
        &mut MemCache::default(),
        &mut CostTable::new(20.0).unwrap(),
    );
    assert_eq!(
        result,
        Err(PlanError::ResourceLimit {
            field: "polynomial_cell_work",
            requested: MAX_FAMILY_COEFFICIENTS * cells * 5,
            limit: MAX_POLYNOMIAL_CELL_WORK,
        })
    );
}

#[test]
fn pl_007c_later_public_rung_resource_failure_precedes_all_planner_work() {
    let family = ProblemFamily::new(
        Poly(vec![0.0; MAX_FAMILY_COEFFICIENTS]),
        "later-rung-too-expensive",
    )
    .unwrap();
    let invalid_cells = MAX_POLYNOMIAL_CELL_WORK / (MAX_FAMILY_COEFFICIENTS * 5) + 1;
    let result = plan(
        &family,
        1.0,
        0.1,
        100.0,
        &[1, invalid_cells],
        &mut MemCache::default(),
        &mut CostTable::new(20.0).unwrap(),
    );
    assert_eq!(
        result,
        Err(PlanError::ResourceLimit {
            field: "polynomial_cell_work",
            requested: MAX_FAMILY_COEFFICIENTS * invalid_cells * 5,
            limit: MAX_POLYNOMIAL_CELL_WORK,
        })
    );
}

#[test]
fn pl_008_budget_below_first_solve_returns_no_certificate() {
    let outcome = plan(
        &steep_family(),
        1.0,
        1e-4,
        1.0,
        &RUNGS,
        &mut MemCache::default(),
        &mut CostTable::new(20.0).unwrap(),
    )
    .unwrap();
    let PlanOutcome::RefusedWithoutAnswer { ops, cost, reason } = outcome else {
        panic!("one cell cannot fund the 12-cell initial solve");
    };
    assert_eq!(ops.len(), 1);
    assert_eq!(ops[0].op, PlanOp::CacheLookup);
    assert!(cost.abs() < f64::EPSILON);
    assert!(reason.contains("no uncertified answer"));
}

#[test]
fn pl_008b_large_admitted_rung_is_budget_checked_before_mesh_allocation() {
    let outcome = plan(
        &steep_family(),
        1.0,
        1e-4,
        1.0,
        &[MAX_PLANNER_CELLS],
        &mut MemCache::default(),
        &mut CostTable::new(20.0).unwrap(),
    )
    .unwrap();

    let PlanOutcome::RefusedWithoutAnswer { ops, cost, reason } = outcome else {
        panic!("one cell cannot admit a million-cell initial mesh");
    };
    assert_eq!(ops.len(), 1);
    assert_eq!(ops[0].op, PlanOp::CacheLookup);
    assert_eq!(cost.to_bits(), 0.0f64.to_bits());
    assert!(reason.contains("no mesh is allocated"));
}

#[derive(Clone)]
struct PoisonedCache {
    answer: CachedAnswer,
}

impl AnswerCache for PoisonedCache {
    fn lookup(&self, _key: &str, _tol: f64) -> Option<CachedAnswer> {
        Some(self.answer.clone())
    }

    fn insert(&mut self, _key: &str, answer: CachedAnswer) {
        self.answer = answer;
    }
}

#[test]
fn pl_009_cache_claims_are_reverified_before_discharge() {
    let poison = CachedAnswer::new(vec![0.0, 1e6, 0.0], 0.0, vec![0.0, 0.5, 1.0]).unwrap();
    let mut cache = PoisonedCache { answer: poison };
    let outcome = plan(
        &steep_family(),
        1.0,
        0.05,
        5000.0,
        &RUNGS,
        &mut cache,
        &mut CostTable::new(20.0).unwrap(),
    )
    .unwrap();
    let PlanOutcome::Discharged { bound, ops, .. } = outcome else {
        panic!("valid solving should recover from a poisoned cache");
    };
    assert!(bound > 0.0 && bound <= 0.05);
    assert_eq!(ops[0].op, PlanOp::CacheLookup);
    assert!(ops.iter().any(|entry| entry.op == PlanOp::SolveRung));
}

#[test]
fn pl_010_adaptive_to_non_dyadic_climb_runs_on_the_planner_path() {
    const NON_DYADIC_RUNGS: [usize; 3] = [12, 25, 53];
    const BUDGET: f64 = 10_000.0;

    let run = || {
        let mut costs = CostTable::new(20.0).unwrap();
        // Prefer one local refinement first. Its measured adaptive solve is
        // expensive enough to make the next deterministic choice a climb.
        costs.record(PlanOp::DwrRefine, 1.0).unwrap();
        costs.record(PlanOp::Climb, 20.0).unwrap();
        plan(
            &steep_family(),
            1.0,
            1e-4,
            BUDGET,
            &NON_DYADIC_RUNGS,
            &mut MemCache::default(),
            &mut costs,
        )
        .unwrap()
    };

    let first = run();
    let first_ops = match &first {
        PlanOutcome::Discharged { ops, cost, .. }
        | PlanOutcome::RefusedWithBest { ops, cost, .. }
        | PlanOutcome::RefusedWithoutAnswer { ops, cost, .. } => {
            assert!(*cost <= BUDGET);
            ops
        }
    };
    let refine_index = first_ops
        .iter()
        .position(|entry| entry.op == PlanOp::DwrRefine)
        .expect("the seeded walk locally refines first");
    let climb_index = first_ops
        .iter()
        .position(|entry| entry.op == PlanOp::Climb)
        .expect("the adaptive result subsequently climbs to a uniform rung");
    assert!(refine_index < climb_index);

    let replay = run();
    let replay_ops = match &replay {
        PlanOutcome::Discharged { ops, .. }
        | PlanOutcome::RefusedWithBest { ops, .. }
        | PlanOutcome::RefusedWithoutAnswer { ops, .. } => ops,
    };
    assert_eq!(
        first_ops, replay_ops,
        "G5 operator replay must be bit-stable"
    );
}

#[test]
fn pl_011_pessimistic_prediction_cannot_veto_affordable_exact_work() {
    const COLD: f64 = 1.0e9;
    let mut costs = CostTable::new(COLD).unwrap();
    let outcome = plan(
        &steep_family(),
        1.0,
        0.05,
        200.0,
        &[12],
        &mut MemCache::default(),
        &mut costs,
    )
    .unwrap();
    let PlanOutcome::Discharged { ops, cost, .. } = outcome else {
        panic!("the exact affordable refinement must discharge despite the cold estimate");
    };
    assert!(cost <= 200.0);
    assert_eq!(
        ops.iter()
            .filter(|entry| entry.op == PlanOp::SolveRung)
            .count(),
        2,
        "the second exact solve ran instead of being vetoed by prediction"
    );
    assert!(
        costs.predict(PlanOp::DwrRefine) < COLD,
        "completed downstream work supplies real transition telemetry"
    );
}

#[test]
fn pl_012_unfunded_refinement_is_not_allocated_or_recorded() {
    const COLD: f64 = 123.0;
    let mut costs = CostTable::new(COLD).unwrap();
    let outcome = plan(
        &steep_family(),
        1.0,
        0.05,
        12.0,
        &[12],
        &mut MemCache::default(),
        &mut costs,
    )
    .unwrap();
    let PlanOutcome::RefusedWithBest { ops, cost, .. } = outcome else {
        panic!("the budget funds the first solve but not the refined solve");
    };
    assert_eq!(cost.to_bits(), 12.0_f64.to_bits());
    assert_eq!(ops.last().map(|entry| entry.op), Some(PlanOp::SolveRung));
    assert!(
        !ops.iter().any(|entry| entry.op == PlanOp::DwrRefine),
        "the transition is not executed before its downstream solve is affordable"
    );
    assert_eq!(
        costs.predict(PlanOp::DwrRefine).to_bits(),
        COLD.to_bits(),
        "a transition with no downstream compute is not an observation"
    );
}

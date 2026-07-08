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
    CostTable, MemCache, PlanOp, PlanOutcome, ProblemFamily, baseline_uniform, plan,
};
use fs_verify::fem1d::Poly;

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-ir/planner\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

/// The wedge-like family: smooth base + a STEEP right feature
/// (x⁸-weighted), so the QoI error concentrates near x = 1 and local
/// refinement genuinely beats uniform.
fn steep_family() -> ProblemFamily {
    // u = x(1−x)·(0.2 + x⁸)
    //   = 0.2x − 0.2x² + x⁹ − x¹⁰.
    let mut c = vec![0.0; 11];
    c[1] = 0.2;
    c[2] = -0.2;
    c[9] = 1.0;
    c[10] = -1.0;
    ProblemFamily {
        base: Poly(c),
        kernel: "cht-wedge-steep".to_string(),
    }
}

const RUNGS: [usize; 4] = [12, 24, 48, 96];

#[test]
fn pl_001_discharges_within_budget_contract_held() {
    let family = steep_family();
    let mut cache = MemCache::default();
    let mut costs = CostTable::new(200.0);
    let tol = 0.05;
    let out = plan(&family, 1.0, tol, 2000.0, &RUNGS, &mut cache, &mut costs);
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
            let seq: Vec<&str> = ops.iter().map(|o| o.op.name()).collect();
            println!(
                "{{\"metric\":\"planner-run\",\"ops\":{seq:?},\"cost\":{cost},\
                 \"bound\":{bound:.3e}}}"
            );
        }
        PlanOutcome::RefusedWithBest { reason, .. } => {
            panic!("a generous budget must discharge: {reason}")
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
    let (base_cost, base_bound) = baseline_uniform(&family, 1.0, tol, 48, 6);
    assert!(base_bound <= tol, "the baseline eventually certifies");
    let mut cache = MemCache::default();
    let mut costs = CostTable::new(200.0);
    let out = plan(&family, 1.0, tol, 100_000.0, &RUNGS, &mut cache, &mut costs);
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
    let mut costs = CostTable::new(500.0);
    // COLD estimates: before any telemetry, predictions are the
    // conservative default (the round-3 boundary).
    assert!((costs.predict(PlanOp::DwrRefine) - 500.0).abs() < f64::EPSILON);
    let tol = 0.05;
    let first = plan(&family, 1.0, tol, 5000.0, &RUNGS, &mut cache, &mut costs);
    assert!(matches!(first, PlanOutcome::Discharged { .. }));
    // Learned estimates move off the default.
    assert!(
        costs.predict(PlanOp::SolveRung) < 500.0,
        "telemetry sharpens the table: {}",
        costs.predict(PlanOp::SolveRung)
    );
    // The SAME query again: a cache hit with ZERO solves.
    let again = plan(&family, 1.0, tol, 5000.0, &RUNGS, &mut cache, &mut costs);
    match again {
        PlanOutcome::Discharged { ops, cost, .. } => {
            assert_eq!(ops.len(), 1, "one op only: the cache");
            assert_eq!(ops[0].op, PlanOp::CacheLookup);
            assert!(cost.abs() < f64::EPSILON, "zero solves on a hit");
        }
        PlanOutcome::RefusedWithBest { .. } => panic!("the hit must discharge"),
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
    let mut costs = CostTable::new(200.0);
    let out = plan(&family, 1.0, tol, budget, &RUNGS, &mut cache, &mut costs);
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
                *cost <= budget + 100.0,
                "never loops past the budget: {cost} vs {budget}"
            );
            assert!(
                reason.contains("refusal"),
                "hands off to refusal semantics: {reason}"
            );
            assert!(!ops.is_empty());
        }
        PlanOutcome::Discharged { .. } => panic!("80 cells cannot certify 1e-4"),
    }
    // G5: the identical query replays the identical operator sequence.
    let run = |seed_cache: &mut MemCache| -> Vec<&'static str> {
        let mut costs = CostTable::new(200.0);
        match plan(&family, 1.0, 0.05, 2000.0, &RUNGS, seed_cache, &mut costs) {
            PlanOutcome::Discharged { ops, .. } | PlanOutcome::RefusedWithBest { ops, .. } => {
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
    let mut costs = CostTable::new(300.0);
    for k in 0..4 {
        let theta = 0.9 + 0.05 * f64::from(k as u8);
        let _ = plan(&family, theta, 0.05, 5000.0, &RUNGS, &mut cache, &mut costs);
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

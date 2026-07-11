//! Anytime + refusal conformance (the lmp4.17 bead; runs under the
//! `ladder-planner` feature). Acceptance: every query returns an
//! immediate colored interval that tightens MONOTONICALLY with budget;
//! each result carries a valid priced "what would tighten this" hint;
//! an under-budget query is REFUSED with the achieved interval and the
//! price of the gap — never a silent point estimate; replays reproduce
//! the trajectory (G5).
#![cfg(feature = "ladder-planner")]

use std::cell::Cell;
use std::rc::Rc;

use fs_evidence::Color;
use fs_ir::anytime::{AnytimeControl, MAX_BUDGET_RUNGS, run_anytime, run_anytime_observed};
use fs_ir::planner::{
    AnswerCache, CachedAnswer, CostTable, MemCache, PlanError, PlanOp, ProblemFamily,
};
use fs_verify::fem1d::Poly;

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-ir/anytime\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn steep_family() -> ProblemFamily {
    let mut c = vec![0.0; 6];
    c[1] = 0.2;
    c[2] = -0.2;
    c[4] = 1.0;
    c[5] = -1.0;
    ProblemFamily::new(Poly(c), "cht-wedge-steep").unwrap()
}

const RUNGS: [usize; 4] = [12, 24, 48, 96];

#[test]
fn an_001_immediate_interval_monotone_tightening() {
    let family = steep_family();
    let mut cache = MemCache::default();
    let mut costs = CostTable::new(200.0).unwrap();
    let budgets = [15.0, 40.0, 120.0, 400.0];
    let report = run_anytime(&family, 1.0, 5e-3, &budgets, &RUNGS, &mut cache, &mut costs).unwrap();
    // Immediate: the FIRST rung already returned a certified interval.
    assert!(!report.trajectory.is_empty());
    let first = &report.trajectory[0];
    assert!(first.bound.is_finite(), "an immediate wide interval exists");
    assert!(
        matches!(first.color, Color::Verified { .. }),
        "the operator knows what kind of answer they hold: {:?}",
        first.color
    );
    // Monotone tightening with budget.
    for pair in report.trajectory.windows(2) {
        assert!(
            pair[1].bound <= pair[0].bound + 1e-15,
            "intervals tighten monotonically: {} -> {}",
            pair[0].bound,
            pair[1].bound
        );
    }
    // G5: an identical replay reproduces the identical trajectory.
    let replay = run_anytime(
        &family,
        1.0,
        5e-3,
        &budgets,
        &RUNGS,
        &mut MemCache::default(),
        &mut CostTable::new(200.0).unwrap(),
    )
    .unwrap();
    assert_eq!(replay.trajectory.len(), report.trajectory.len());
    for (a, b) in report.trajectory.iter().zip(&replay.trajectory) {
        assert_eq!(a.bound.to_bits(), b.bound.to_bits(), "bit-equal trajectory");
        assert_eq!(a.hint, b.hint);
        assert_eq!(a.verifier_family, b.verifier_family);
        assert_eq!(a.flux_hash, b.flux_hash);
        assert_eq!(a.color, b.color);
    }
    println!(
        "{{\"metric\":\"anytime-trajectory\",\"bounds\":{:?}}}",
        report
            .trajectory
            .iter()
            .map(|s| s.bound)
            .collect::<Vec<_>>()
    );
    verdict(
        "an-001",
        "immediate verified interval; monotone tightening across the budget ladder; \
         bit-equal replay (G5)",
    );
}

#[test]
fn an_002_refusal_teaches_with_the_gap_price() {
    let family = steep_family();
    let mut cache = MemCache::default();
    let mut costs = CostTable::new(200.0).unwrap();
    // A tolerance the small budget ladder cannot reach.
    let report = run_anytime(
        &family,
        1.0,
        1e-4,
        &[20.0, 45.0, 80.0],
        &RUNGS,
        &mut cache,
        &mut costs,
    )
    .unwrap();
    assert!(!report.discharged(), "the query could not discharge");
    let refusal = report.refusal.as_ref().expect("a refusal note exists");
    // The refusal carries: the ACHIEVED interval…
    assert!(refusal.contains("achieved a certified"), "{refusal}");
    assert!(refusal.contains("±"), "the interval is stated: {refusal}");
    // …the PRICE of the gap…
    assert!(
        refusal.contains("more cells"),
        "the gap is priced: {refusal}"
    );
    // …and NO silent point estimate.
    assert!(
        refusal.contains("No best-effort point estimate"),
        "the honesty clause is explicit: {refusal}"
    );
    // Every trajectory step still carried a certified interval + color.
    for step in &report.trajectory {
        assert!(step.bound.is_finite());
        assert!(matches!(step.color, Color::Verified { .. }));
        assert!(!step.discharged);
    }
    println!("{{\"metric\":\"refusal\",\"note\":{refusal:?}}}");
    verdict(
        "an-002",
        "the impossible query refuses with the achieved interval, the priced gap, and \
         the explicit no-point-estimate clause",
    );
}

#[test]
fn an_003_hint_names_a_real_move_and_the_hot_region() {
    let family = steep_family();
    let mut cache = MemCache::default();
    let mut costs = CostTable::new(200.0).unwrap();
    let report = run_anytime(
        &family,
        1.0,
        1e-4,
        &[60.0, 100.0],
        &RUNGS,
        &mut cache,
        &mut costs,
    )
    .unwrap();
    let step = report.trajectory.last().expect("steps exist");
    // The hint names a REAL operator from the menu…
    assert!(
        step.hint.contains("dwr-refine") || step.hint.contains("climb"),
        "a real menu move: {}",
        step.hint
    );
    // …prices it…
    assert!(step.hint.contains("more cells"), "priced: {}", step.hint);
    // …and, after local refinement, names WHERE the money goes (the
    // steep feature lives near x = 1).
    if step.hint.contains("region") {
        assert!(
            step.hint.contains("x ∈ ["),
            "the hot region is an interval: {}",
            step.hint
        );
    }
    // Cold-telemetry degradation: a fresh table still yields a priced
    // hint (the generic form).
    let cold_hint =
        fs_ir::anytime::tighten_hint(0.05, 1e-3, 30.0, &CostTable::new(500.0).unwrap(), None)
            .unwrap();
    assert!(
        cold_hint.contains("more cells"),
        "cold hint priced: {cold_hint}"
    );
    verdict(
        "an-003",
        "hints name a real menu move with a price; local refinement names the hot \
         region; cold telemetry degrades to the generic priced form",
    );
}

#[test]
fn an_004_discharge_ends_the_trajectory_and_caches() {
    let family = steep_family();
    let mut cache = MemCache::default();
    let mut costs = CostTable::new(200.0).unwrap();
    let budgets = [15.0, 400.0, 4000.0];
    let report = run_anytime(&family, 1.0, 5e-3, &budgets, &RUNGS, &mut cache, &mut costs).unwrap();
    assert!(report.discharged());
    let last = report.trajectory.last().expect("steps");
    assert!(last.discharged);
    assert!(
        report.trajectory.len() < budgets.len() || last.discharged,
        "the trajectory stops at discharge (no wasted rungs after success)"
    );
    assert!(last.hint.contains("spend nothing"), "{}", last.hint);
    // The follow-up identical query is a pure cache hit at ANY budget.
    let again = run_anytime(&family, 1.0, 5e-3, &[1.0], &RUNGS, &mut cache, &mut costs).unwrap();
    assert!(again.discharged(), "the cached answer discharges instantly");
    assert_eq!(again.trajectory.len(), 1);
    verdict(
        "an-004",
        "discharge terminates the ladder with the 'spend nothing' hint; the repeat \
         query discharges from cache at a 1-cell budget",
    );
}

#[test]
fn an_005_invalid_ladders_and_scalars_refuse_before_planning() {
    let family = steep_family();
    let run = |theta, tol, budgets: &[f64], rungs: &[usize]| {
        run_anytime(
            &family,
            theta,
            tol,
            budgets,
            rungs,
            &mut MemCache::default(),
            &mut CostTable::new(20.0).unwrap(),
        )
    };
    assert_eq!(
        run(1.0, 0.1, &[], &RUNGS),
        Err(PlanError::EmptySequence {
            field: "budget_ladder",
        })
    );
    assert!(matches!(
        run(1.0, 0.1, &[0.0], &RUNGS),
        Err(PlanError::InvalidSequenceEntry {
            field: "budget_ladder",
            index: 0,
            ..
        })
    ));
    assert!(matches!(
        run(1.0, 0.1, &[10.0, 5.0], &RUNGS),
        Err(PlanError::InvalidSequenceEntry {
            field: "budget_ladder",
            index: 1,
            ..
        })
    ));
    assert!(matches!(
        run(1.0, 0.1, &[10.0], &[]),
        Err(PlanError::EmptySequence {
            field: "rung_cells"
        })
    ));
    let final_oversized_budget = u32::try_from(MAX_BUDGET_RUNGS + 1).expect("test cap fits u32");
    let oversized_budgets = (1..=final_oversized_budget)
        .map(f64::from)
        .collect::<Vec<_>>();
    assert_eq!(
        run(1.0, 0.1, &oversized_budgets, &RUNGS),
        Err(PlanError::ResourceLimit {
            field: "budget_ladder",
            requested: MAX_BUDGET_RUNGS + 1,
            limit: MAX_BUDGET_RUNGS,
        })
    );
    assert!(matches!(
        run(f64::NAN, 0.1, &[10.0], &RUNGS),
        Err(PlanError::InvalidScalar { field: "theta", .. })
    ));
    assert!(matches!(
        run(1.0, 0.0, &[10.0], &RUNGS),
        Err(PlanError::InvalidScalar {
            field: "tolerance",
            ..
        })
    ));
}

#[test]
fn an_006_unaffordable_first_solve_carries_no_fake_interval_or_color() {
    let report = run_anytime(
        &steep_family(),
        1.0,
        1e-4,
        &[1.0],
        &RUNGS,
        &mut MemCache::default(),
        &mut CostTable::new(20.0).unwrap(),
    )
    .unwrap();
    assert!(report.trajectory.is_empty());
    assert!(!report.discharged());
    let refusal = report.refusal.expect("explicit no-certificate refusal");
    assert!(refusal.contains("without a certified interval"));
    assert!(refusal.contains("No best-effort point estimate"));

    let costs = CostTable::new(20.0).unwrap();
    assert!(fs_ir::anytime::tighten_hint(f64::NAN, 0.1, 1.0, &costs, None).is_err());
    assert!(fs_ir::anytime::tighten_hint(1.0, 0.0, 1.0, &costs, None).is_err());
    assert!(fs_ir::anytime::tighten_hint(1.0, 0.1, -1.0, &costs, None).is_err());
    assert!(fs_ir::anytime::tighten_hint(1.0, 0.1, 1.0, &costs, Some((0.8, 0.2))).is_err());
}

#[derive(Default)]
struct CountingCache {
    inner: MemCache,
    lookups: Cell<usize>,
}

impl AnswerCache for CountingCache {
    fn lookup(&self, key: &str, tol: f64) -> Option<CachedAnswer> {
        self.lookups.set(self.lookups.get() + 1);
        self.inner.lookup(key, tol)
    }

    fn insert(&mut self, key: &str, answer: CachedAnswer) {
        self.inner.insert(key, answer);
    }
}

#[test]
fn an_007_budget_rungs_replay_one_execution_prefix() {
    let budgets = [1.0, 11.0, 12.0, 13.0, 200.0];
    let mut cache = CountingCache::default();
    let report = run_anytime(
        &steep_family(),
        1.0,
        0.05,
        &budgets,
        &RUNGS,
        &mut cache,
        &mut CostTable::new(1.0e9).unwrap(),
    )
    .unwrap();

    assert_eq!(cache.lookups.get(), 1, "the planner executes only once");
    assert_eq!(
        report.trajectory.first().map(|step| step.budget),
        Some(12.0),
        "no colored point appears before the first 12-cell solve is affordable"
    );
    assert_eq!(
        report.trajectory[0].bound.to_bits(),
        report.trajectory[1].bound.to_bits()
    );
    for pair in report.trajectory.windows(2) {
        assert!(
            pair[1].bound <= pair[0].bound,
            "a single execution prefix can only retain or improve its certificate"
        );
    }
    assert!(report.discharged());
}

#[derive(Default)]
struct SharedCacheCounts {
    lookups: Cell<usize>,
    inserts: Cell<usize>,
}

struct ObservedCache {
    counts: Rc<SharedCacheCounts>,
}

impl AnswerCache for ObservedCache {
    fn lookup(&self, _key: &str, _tol: f64) -> Option<CachedAnswer> {
        self.counts.lookups.set(self.counts.lookups.get() + 1);
        None
    }

    fn insert(&mut self, _key: &str, _answer: CachedAnswer) {
        self.counts.inserts.set(self.counts.inserts.get() + 1);
    }
}

#[test]
fn an_008_observer_stop_prevents_later_work_telemetry_and_cache_insertion() {
    const COLD: f64 = 1.0e9;
    let counts = Rc::new(SharedCacheCounts::default());
    let mut cache = ObservedCache {
        counts: Rc::clone(&counts),
    };
    let mut costs = CostTable::new(COLD).unwrap();
    let mut callbacks = 0usize;
    let report = run_anytime_observed(
        &steep_family(),
        1.0,
        1e-4,
        &[12.0, 400.0],
        &RUNGS,
        &mut cache,
        &mut costs,
        |_| {
            callbacks += 1;
            AnytimeControl::Stop
        },
    )
    .unwrap();

    assert_eq!(callbacks, 1);
    assert!(report.stopped_by_observer());
    assert_eq!(report.trajectory.len(), 1);
    assert_eq!(report.trajectory[0].budget.to_bits(), 12.0_f64.to_bits());
    assert_eq!(report.trajectory[0].spent.to_bits(), 12.0_f64.to_bits());
    assert_eq!(counts.lookups.get(), 1);
    assert_eq!(counts.inserts.get(), 0);
    assert_eq!(costs.predict(PlanOp::DwrRefine).to_bits(), COLD.to_bits());
    assert_eq!(costs.predict(PlanOp::Climb).to_bits(), COLD.to_bits());
    assert_eq!(costs.predict(PlanOp::Speculate).to_bits(), COLD.to_bits());
    let refusal = report.refusal.expect("observer stop is explicit");
    assert!(refusal.contains("12.0 cells spent at budget rung 12.0"));

    let discharge_counts = Rc::new(SharedCacheCounts::default());
    let mut discharge_cache = ObservedCache {
        counts: Rc::clone(&discharge_counts),
    };
    let discharged = run_anytime_observed(
        &steep_family(),
        1.0,
        1.0,
        &[12.0],
        &RUNGS,
        &mut discharge_cache,
        &mut CostTable::new(COLD).unwrap(),
        |_| AnytimeControl::Stop,
    )
    .unwrap();
    assert!(discharged.discharged());
    assert!(discharged.stopped_by_observer());
    assert_eq!(
        discharge_counts.inserts.get(),
        0,
        "the callback observes discharge before the later cache insertion"
    );
}

#[test]
fn an_009_callbacks_are_prefix_ordered_and_hints_exclude_future_telemetry() {
    let mut callback_prefix = Vec::new();
    let mut costs = CostTable::new(1.0).unwrap();
    let report = run_anytime_observed(
        &steep_family(),
        1.0,
        0.05,
        &[12.0, 20.0, 100.0, 400.0],
        &RUNGS,
        &mut MemCache::default(),
        &mut costs,
        |step| {
            callback_prefix.push((
                step.budget.to_bits(),
                step.spent.to_bits(),
                step.hint.clone(),
            ));
            AnytimeControl::Continue
        },
    )
    .unwrap();

    assert_eq!(callback_prefix.len(), report.trajectory.len());
    for (observed, retained) in callback_prefix.iter().zip(&report.trajectory) {
        assert_eq!(observed.0, retained.budget.to_bits());
        assert_eq!(observed.1, retained.spent.to_bits());
        assert_eq!(observed.2, retained.hint);
    }
    assert!(
        report.trajectory[0].hint.contains("dwr-refine"),
        "the first hint sees the contemporaneous cold tie, not later refinement telemetry: {}",
        report.trajectory[0].hint
    );
    assert_eq!(report.trajectory[0].verifier_family, "equilibrated-flux-1d");
    assert_eq!(
        report.trajectory[0].color,
        Color::Verified {
            lo: 0.0,
            hi: report.trajectory[0].bound,
        }
    );
    for pair in report.trajectory.windows(2) {
        assert!(pair[0].budget < pair[1].budget);
        assert!(pair[0].spent <= pair[1].spent);
        assert!(pair[0].bound >= pair[1].bound);
    }
}

#[test]
fn an_010_speculation_stop_retains_current_but_not_later_telemetry() {
    let mut costs = CostTable::new(50.0).unwrap();
    costs.record(PlanOp::DwrRefine, 1_000.0).unwrap();
    costs.record(PlanOp::Climb, 1.0).unwrap();
    let mut callbacks = 0usize;
    let report = run_anytime_observed(
        &steep_family(),
        1.0,
        1e-4,
        &[12.0, 17.0],
        &[12, 24],
        &mut MemCache::default(),
        &mut costs,
        |_| {
            callbacks += 1;
            if callbacks == 2 {
                AnytimeControl::Stop
            } else {
                AnytimeControl::Continue
            }
        },
    )
    .unwrap();

    assert_eq!(callbacks, 2);
    assert!(report.stopped_by_observer());
    assert_eq!(report.trajectory[1].budget.to_bits(), 17.0_f64.to_bits());
    assert!(report.trajectory[1].spent < 17.0);
    assert!(
        costs.predict(PlanOp::Climb) < 3.0,
        "completed speculative transition is visible at its certificate checkpoint"
    );
    assert!((costs.predict(PlanOp::Speculate) - 4.8).abs() < 1e-12);
    assert_eq!(
        costs.predict(PlanOp::SolveRung).to_bits(),
        12.0_f64.to_bits(),
        "the later fine-rung solve never ran or entered telemetry"
    );
}

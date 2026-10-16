//! fs-plan moonshot battery (bead gp3.9 V2, feature
//! `moonshot-planner`): the co-optimizer is EXACT against the brute
//! oracle, the convex water-filling agrees with CMA-ES where the
//! model is convex, CMA-ES escapes where it is not, and the
//! fixture-matrix scoreboard (V2 vs V1 vs hand) is generated — the
//! promotion-gate evidence SHAPE, with the flagship-set gate honestly
//! left to the huq.15 Gauntlet.

use fs_plan::moonshot::{
    RateKnob, ScoreRow, cma_continuous, optimize_exact, rate_error, waterfill,
};
use fs_plan::{AllocProblem, AllocationError, Knob, KnobSetting, allocate, oracle_min_error};

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

fn ks(label: &str, error: f64, cost: f64) -> KnobSetting {
    KnobSetting {
        label: label.to_owned(),
        error,
        cost,
    }
}

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / (1u64 << 53) as f64
}

fn random_problem(seed: &mut u64, tracks: usize) -> AllocProblem {
    let mut knobs = Vec::new();
    for k in 0..4 {
        let mut settings = Vec::new();
        let mut err = 1.0 + lcg(seed);
        let mut cost = 0.5 + lcg(seed);
        let nset = 2 + (lcg(seed) * 3.0) as usize;
        for s in 0..nset {
            settings.push(ks(&format!("k{k}s{s}"), err, cost));
            err *= 0.25 + 0.3 * lcg(seed);
            cost *= 1.8 + 1.5 * lcg(seed);
        }
        knobs
            .push(Knob::new(&format!("knob{k}"), k % tracks, settings).expect("valid random knob"));
    }
    AllocProblem {
        knobs,
        budget_s: 6.0,
        error_target: 0.0,
    }
}

/// ms-001: the DP co-optimizer matches the brute-force oracle exactly
/// (single- and multi-track fixtures) and never loses to greedy V1.
#[test]
fn ms_001_dp_is_exact() {
    let mut seed = 0x300D_u64;
    let mut v2_wins_or_ties = 0usize;
    let mut cases = 0usize;
    for tracks in [1usize, 2] {
        for _ in 0..40 {
            let p = random_problem(&mut seed, tracks);
            let (_, oracle_err) = oracle_min_error(&p)
                .expect("valid fixture")
                .expect("cheapest fits");
            let v2 = optimize_exact(&p, 4000).expect("DP finds a plan");
            // Grid-conservative: DP may only lose oracle by bucket
            // rounding, never win.
            assert!(
                v2.total_error + 1e-12 >= oracle_err,
                "DP beat the oracle?! {} < {oracle_err}",
                v2.total_error
            );
            let rounding = (v2.total_error - oracle_err) / oracle_err.max(1e-30);
            assert!(
                rounding < 0.02,
                "DP bucket rounding too coarse: {rounding:.3} at tracks={tracks}"
            );
            let v1_err = match allocate(&p) {
                Ok(plan) => plan.total_error,
                Err(AllocationError::BudgetInfeasible(inf)) => inf.best_error_in_budget,
                Err(other) => panic!("valid fixture produced unexpected refusal: {other}"),
            };
            if v2.total_error <= v1_err + 1e-12 {
                v2_wins_or_ties += 1;
            }
            cases += 1;
        }
    }
    verdict(
        "ms-001-dp-exact",
        v2_wins_or_ties == cases,
        &format!(
            "{cases} fixtures (1 and 2 tracks): DP within 2% of oracle (bucket rounding only) and never loses to greedy"
        ),
    );
}

/// ms-002: convex cross-check — CMA-ES on the continuous rate-based
/// allocation agrees with the water-filling KKT solution.
#[test]
fn ms_002_cma_agrees_with_waterfill_on_convex() {
    let knobs = [
        RateKnob { a: 0.9, p: 0.7 },
        RateKnob { a: 0.5, p: 1.1 },
        RateKnob { a: 0.2, p: 1.4 },
        RateKnob { a: 0.3, p: 0.9 },
    ];
    let budget = 100.0;
    let wf = waterfill(&knobs, budget);
    let e_wf = rate_error(&knobs, &wf);
    let mut model = |alloc: &[f64]| rate_error(&knobs, alloc);
    let cma = cma_continuous(knobs.len(), budget, &mut model, 0x5EED);
    let e_cma = rate_error(&knobs, &cma);
    let rel = (e_cma - e_wf).abs() / e_wf;
    verdict(
        "ms-002-cma-vs-waterfill",
        rel < 5e-3 && e_cma + 1e-12 >= e_wf * 0.999,
        &format!("convex model: waterfill {e_wf:.6}, CMA-ES {e_cma:.6} (rel gap {rel:.2e})"),
    );
}

/// ms-003: NON-convex model — a surrogate with an activation discount
/// (error drops only past a spend threshold). Water-filling's smooth
/// KKT never sees the cliff; CMA-ES finds it.
#[test]
fn ms_003_cma_escapes_nonconvex() {
    let smooth = [RateKnob { a: 0.8, p: 0.8 }, RateKnob { a: 0.6, p: 0.9 }];
    let budget = 60.0;
    // Knob 2 (index 2): surrogate — error 0.5 below 30s spend, 0.02
    // at/above (training threshold): a cliff, not a rate curve.
    let mut model = |alloc: &[f64]| -> f64 {
        let base = rate_error(&smooth, &alloc[..2]);
        let surr = if alloc[2] >= 30.0 { 0.02 } else { 0.5 };
        base + surr
    };
    // Water-filling can only handle the smooth part: model the cliff
    // as its smooth secant (the natural convexification) and round.
    let wf_smooth = waterfill(
        &[
            smooth[0],
            smooth[1],
            RateKnob {
                a: 0.5 * 30.0,
                p: 1.0,
            }, // secant through the cliff
        ],
        budget,
    );
    let e_wf = model(&wf_smooth);
    let cma = cma_continuous(3, budget, &mut model, 0xC11FF);
    let e_cma = model(&cma);
    verdict(
        "ms-003-cma-nonconvex",
        e_cma < e_wf * 0.8 && cma[2] >= 30.0,
        &format!(
            "cliff model: convexified waterfill {e_wf:.4}, CMA-ES {e_cma:.4} (surrogate spend {:.1}s >= 30s threshold)",
            cma[2]
        ),
    );
}

/// ms-004: the scoreboard — V2 vs V1 vs HAND (uniform mid-settings)
/// across the fixture matrix; accuracy-per-second rows ledgered. The
/// flagship promotion gate remains huq.15's call (feature ships OFF).
#[test]
fn ms_004_scoreboard_v2_beats_hand_on_fixtures() {
    let mut seed = 0xB0A2D_u64;
    let mut rows = Vec::new();
    for i in 0..25 {
        let p = random_problem(&mut seed, 1 + i % 2);
        let hand_choice: Vec<usize> = p.knobs.iter().map(|k| k.settings.len() / 2).collect();
        let hand_wall =
            fs_plan::alloc::plan_wall_clock(&p.knobs, &hand_choice).expect("valid hand choice");
        let hand_error = if hand_wall <= p.budget_s {
            fs_plan::alloc::plan_total_error(&p.knobs, &hand_choice).expect("valid hand choice")
        } else {
            // Hand plan busts the budget: charge it the cheapest plan
            // (the generous reading).
            fs_plan::alloc::plan_total_error(&p.knobs, &vec![0; p.knobs.len()])
                .expect("valid cheapest choice")
        };
        let v1_error = match allocate(&p) {
            Ok(plan) => plan.total_error,
            Err(AllocationError::BudgetInfeasible(inf)) => inf.best_error_in_budget,
            Err(other) => panic!("valid fixture produced unexpected refusal: {other}"),
        };
        let v2_error = optimize_exact(&p, 4000).expect("plan").total_error;
        rows.push(ScoreRow {
            fixture: format!("fixture-{i}"),
            hand_error,
            v1_error,
            v2_error,
        });
    }
    let wins = rows.iter().filter(|r| r.v2_wins()).count();
    for r in &rows {
        println!(
            "{{\"scoreboard\":\"{}\",\"hand\":{:.4},\"v1\":{:.4},\"v2\":{:.4}}}",
            r.fixture, r.hand_error, r.v1_error, r.v2_error
        );
    }
    verdict(
        "ms-004-scoreboard",
        wins == rows.len(),
        &format!(
            "V2 beats-or-ties V1 AND hand allocation on {wins}/{} fixtures (flagship promotion stays with huq.15 — feature ships OFF)",
            rows.len()
        ),
    );
}

/// ms-006: a setting whose cost dwarfs the bucket step must not overflow the
/// knapsack DP — it simply never fits. Regression for `b + buckets` where
/// `buckets` saturated to usize::MAX (debug/test panic; release wrap poisoning
/// a valid bucket with the huge-cost setting → spurious infeasible).
#[test]
fn ms_006_huge_cost_setting_does_not_overflow_the_dp() {
    let knobs = vec![
        Knob::new("a", 0, vec![ks("a0", 1.0, 1.0), ks("a1", 0.5, 2.0)]).expect("knob a"),
        Knob::new("b", 0, vec![ks("b0", 1.0, 1.0), ks("b_huge", 0.1, 1e20)]).expect("knob b"),
    ];
    let problem = AllocProblem {
        knobs,
        budget_s: 6.0,
        error_target: 0.0,
    };
    // Old code panicked here on the usize overflow; the unaffordable setting is
    // now skipped and a plan over the affordable settings is returned.
    let plan = optimize_exact(&problem, 64);
    verdict(
        "ms-006",
        plan.is_some(),
        "1e20-cost knob is skipped; a feasible plan over affordable settings is returned",
    );
}

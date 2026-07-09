//! fs-plan allocator battery (bead gp3.9 V1): greedy vs the exact
//! oracle, structured infeasibility with verified relaxations, online
//! re-planning that reacts when and only when estimates warrant,
//! tropical slack upgrades, the §11.4 "drag to 2% in 2h" scenario, and
//! determinism.

use fs_plan::{AllocProblem, Allocator, Knob, KnobSetting, allocate, oracle_min_error};

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

/// Random fixture: 4 knobs × up to 4 settings, geometric error decay,
/// random costs, single track.
fn random_problem(seed: &mut u64) -> AllocProblem {
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
        knobs.push(Knob::new(&format!("knob{k}"), 0, settings));
    }
    AllocProblem {
        knobs,
        budget_s: 6.0,
        error_target: 0.0, // "as low as possible" — drives full greedy
    }
}

/// pl-001: greedy is feasible and near-oracle on the random fixture
/// matrix (the exact oracle is brute force at fixture scale).
#[test]
fn pl_001_greedy_near_oracle() {
    let mut seed = 0x9_1e57_u64;
    let mut worst_ratio = 1.0f64;
    for _ in 0..60 {
        let p = random_problem(&mut seed);
        // With target 0 the greedy runs until budget-stalled, then
        // reports infeasible carrying its best-in-budget error.
        let greedy_err = match allocate(&p) {
            Ok(plan) => plan.total_error,
            Err(inf) => inf.best_error_in_budget,
        };
        let (_, oracle_err) = oracle_min_error(&p).expect("cheapest plan always fits");
        assert!(
            greedy_err + 1e-12 >= oracle_err,
            "greedy beats the exact oracle?! {greedy_err} < {oracle_err}"
        );
        worst_ratio = worst_ratio.max(greedy_err / oracle_err.max(1e-30));
    }
    verdict(
        "pl-001-greedy-near-oracle",
        worst_ratio < 2.5,
        &format!("60 random fixtures: worst greedy/oracle error ratio {worst_ratio:.3}"),
    );
}

/// pl-002: structured infeasibility — the relaxations are RANKED and
/// VERIFIED (re-planning at the suggested budget succeeds).
#[test]
fn pl_002_infeasibility_is_structured_and_actionable() {
    let knobs = vec![
        Knob::new(
            "mesh",
            0,
            vec![ks("coarse", 1.0, 1.0), ks("fine", 0.1, 10.0)],
        ),
        Knob::new("order", 0, vec![ks("p1", 0.5, 1.0), ks("p3", 0.05, 8.0)]),
    ];
    let p = AllocProblem {
        knobs,
        budget_s: 3.0,
        error_target: 0.2,
    };
    let inf = allocate(&p).expect_err("3s cannot reach 0.2");
    let feasible_again = allocate(&AllocProblem {
        budget_s: inf.budget_needed_for_target * 1.0001,
        ..p.clone()
    });
    verdict(
        "pl-002-structured-infeasibility",
        !inf.relaxations.is_empty()
            && inf.best_error_in_budget > p.error_target
            && feasible_again.is_ok(),
        &format!(
            "best-in-budget {:.3}, budget-for-target {:.1}s, {} ranked relaxations; replan at suggested budget: {}",
            inf.best_error_in_budget,
            inf.budget_needed_for_target,
            inf.relaxations.len(),
            feasible_again.is_ok()
        ),
    );
}

/// pl-003: online re-planning — the plan changes when and only when
/// the a-posteriori estimate warrants it.
#[test]
fn pl_003_online_replanning_reacts_correctly() {
    let knobs = vec![
        Knob::new(
            "mesh",
            0,
            vec![
                ks("m16", 0.8, 1.0),
                ks("m32", 0.2, 4.0),
                ks("m64", 0.05, 16.0),
            ],
        ),
        Knob::new(
            "samples",
            0,
            vec![ks("n100", 0.6, 1.0), ks("n1k", 0.06, 6.0)],
        ),
    ];
    let p = AllocProblem {
        knobs,
        budget_s: 12.0,
        error_target: 0.30,
    };
    let mut alloc = Allocator::new(p);
    let plan0 = alloc.replan().expect("baseline feasible");
    // Unwarranted update: a tiny refinement of an estimate that does
    // not change any greedy comparison.
    alloc.observe_error(0, 0, 0.79);
    let plan1 = alloc.replan().expect("still feasible");
    let unchanged = plan0.choice == plan1.choice;
    // Warranted update: the DWR estimate reveals the coarse mesh is
    // far BETTER than modeled — the expensive mesh upgrade is no
    // longer worth it and the plan must change.
    alloc.observe_error(0, 0, 0.10);
    let plan2 = alloc.replan().expect("feasible after update");
    let changed = plan2.choice != plan0.choice;
    verdict(
        "pl-003-online-replanning",
        unchanged && changed,
        &format!(
            "unwarranted update: plan {:?} unchanged = {unchanged}; warranted update: plan {:?} -> {:?}, changed = {changed}",
            plan0.choice, plan0.choice, plan2.choice
        ),
    );
}

/// pl-004: tropical tracks — an upgrade OFF the critical path costs
/// only slack and the planner takes it for free (wall unchanged).
#[test]
fn pl_004_tropical_slack_upgrade() {
    let knobs = vec![
        // Track 0 dominates the wall-clock (the critical path).
        Knob::new("cfd", 0, vec![ks("coarse", 0.5, 10.0)]),
        // Track 1 has slack: upgrading is free in wall-clock.
        Knob::new(
            "render",
            1,
            vec![ks("preview", 0.4, 1.0), ks("final", 0.05, 6.0)],
        ),
    ];
    let p = AllocProblem {
        knobs,
        budget_s: 10.0,
        error_target: 0.6,
    };
    let plan = allocate(&p).expect("feasible");
    verdict(
        "pl-004-tropical-slack",
        plan.choice == vec![0, 1] && (plan.wall_clock - 10.0).abs() < 1e-12,
        &format!(
            "choice {:?}, wall {:.1}s — the render upgrade rode the slack (rationale: {})",
            plan.choice,
            plan.wall_clock,
            plan.rationale.join("; ")
        ),
    );
}

/// pl-005: the §11.4 scenario — "drag to 2% in 2 hours" over the five
/// canonical knobs, rate-based settings; the plan must meet the target
/// inside the budget with a non-empty rationale, deterministically.
#[test]
fn pl_005_drag_to_two_percent_in_two_hours() {
    let ladder = |name: &str, track: usize, a: f64, p_exp: f64, costs: &[f64]| -> Knob {
        let settings = costs
            .iter()
            .enumerate()
            .map(|(i, &c)| ks(&format!("{name}{i}"), a * c.powf(-p_exp), c))
            .collect();
        Knob::new(name, track, settings)
    };
    let knobs = vec![
        ladder("mesh", 0, 0.9, 0.7, &[60.0, 240.0, 960.0, 3840.0]),
        ladder("order", 0, 0.5, 1.1, &[30.0, 120.0, 480.0]),
        ladder("solver-tol", 0, 0.2, 1.4, &[20.0, 80.0, 320.0]),
        ladder("surrogate", 1, 0.3, 0.9, &[100.0, 400.0, 1600.0]),
        ladder("samples", 1, 0.4, 0.5, &[200.0, 800.0, 3200.0]),
    ];
    let p = AllocProblem {
        knobs,
        budget_s: 7200.0,
        error_target: 0.02,
    };
    let plan = allocate(&p).expect("2% in 2h is feasible on this model");
    let plan_b = allocate(&p).expect("determinism replay");
    verdict(
        "pl-005-drag-scenario",
        plan.total_error <= 0.02
            && plan.wall_clock <= 7200.0
            && !plan.rationale.is_empty()
            && plan.choice == plan_b.choice,
        &format!(
            "error {:.4} <= 2% at wall {:.0}s <= 7200s; {} rationale lines; deterministic replay",
            plan.total_error,
            plan.wall_clock,
            plan.rationale.len()
        ),
    );
}

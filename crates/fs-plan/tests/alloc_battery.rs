//! fs-plan allocator battery (bead gp3.9 V1): greedy vs the exact
//! oracle, structured infeasibility with verified relaxations, online
//! re-planning that reacts when and only when estimates warrant,
//! tropical slack upgrades, the §11.4 "drag to 2% in 2h" scenario, and
//! determinism.

use fs_plan::{
    AllocProblem, AllocationError, Allocator, Knob, KnobSetting, MAX_ALLOCATION_KNOBS,
    MAX_EXECUTION_TRACKS, MAX_ORACLE_COMBINATIONS, MAX_SETTINGS_PER_KNOB, MAX_TOTAL_SETTINGS,
    PlanInputError, allocate, oracle_min_error,
};

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
        knobs.push(Knob::new(&format!("knob{k}"), 0, settings).expect("valid random knob"));
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
            Err(AllocationError::BudgetInfeasible(inf)) => inf.best_error_in_budget,
            Err(other) => panic!("valid fixture produced unexpected refusal: {other}"),
        };
        let (_, oracle_err) = oracle_min_error(&p)
            .expect("valid fixture")
            .expect("cheapest plan always fits");
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
        )
        .expect("valid mesh knob"),
        Knob::new("order", 0, vec![ks("p1", 0.5, 1.0), ks("p3", 0.05, 8.0)])
            .expect("valid order knob"),
    ];
    let p = AllocProblem {
        knobs,
        budget_s: 3.0,
        error_target: 0.2,
    };
    let inf = match allocate(&p).expect_err("3s cannot reach 0.2") {
        AllocationError::BudgetInfeasible(inf) => inf,
        other => panic!("expected target infeasibility, got {other}"),
    };
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
        )
        .expect("valid mesh knob"),
        Knob::new(
            "samples",
            0,
            vec![ks("n100", 0.6, 1.0), ks("n1k", 0.06, 6.0)],
        )
        .expect("valid samples knob"),
    ];
    let p = AllocProblem {
        knobs,
        budget_s: 12.0,
        error_target: 0.30,
    };
    let mut alloc = Allocator::new(p).expect("valid allocation problem");
    let plan0 = alloc.replan().expect("baseline feasible");
    // Unwarranted update: a tiny refinement of an estimate that does
    // not change any greedy comparison.
    alloc.observe_error(0, 0, 0.79).expect("finite observation");
    let plan1 = alloc.replan().expect("still feasible");
    let unchanged = plan0.choice == plan1.choice;
    // Warranted update: the DWR estimate reveals the coarse mesh is
    // far BETTER than modeled — the expensive mesh upgrade is no
    // longer worth it and the plan must change.
    alloc.observe_error(0, 0, 0.10).expect("finite observation");
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
        Knob::new("cfd", 0, vec![ks("coarse", 0.5, 10.0)]).expect("valid CFD knob"),
        // Track 1 has slack: upgrading is free in wall-clock.
        Knob::new(
            "render",
            1,
            vec![ks("preview", 0.4, 1.0), ks("final", 0.05, 6.0)],
        )
        .expect("valid render knob"),
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
        Knob::new(name, track, settings).expect("valid rate ladder")
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

/// pl-007: the cheapest complete plan is still subject to the wall
/// budget. Meeting the error target must never bypass feasibility.
#[test]
fn pl_007_success_never_exceeds_wall_budget() {
    let knob =
        Knob::new("single", 0, vec![ks("only", 0.1, 10.0)]).expect("valid single-setting knob");
    for budget_s in [0.0, 1.0, 9.999, 10.0, 11.0] {
        let result = allocate(&AllocProblem {
            knobs: vec![knob.clone()],
            budget_s,
            error_target: 0.2,
        });
        match result {
            Ok(plan) => assert!(
                plan.wall_clock <= budget_s,
                "successful plan used {}s against a {budget_s}s budget",
                plan.wall_clock
            ),
            Err(AllocationError::MinimumPlanExceedsBudget {
                budget_s: refused_budget,
                minimum_wall_s,
            }) => {
                assert_eq!(refused_budget.to_bits(), budget_s.to_bits());
                assert_eq!(minimum_wall_s.to_bits(), 10.0f64.to_bits());
                assert!(budget_s < minimum_wall_s);
            }
            Err(other) => panic!("unexpected refusal for budget {budget_s}: {other}"),
        }
    }
    verdict(
        "pl-007-budget-safe-success",
        true,
        "target satisfaction cannot return an over-budget baseline plan",
    );
}

/// pl-008: every numeric input domain is enforced with typed errors,
/// including finite values whose aggregate would overflow.
#[test]
fn pl_008_invalid_numeric_domains_are_typed_refusals() {
    for (field, error, cost) in [
        ("error", f64::NAN, 1.0),
        ("error", f64::INFINITY, 1.0),
        ("error", -1.0, 1.0),
        ("cost", 1.0, f64::NAN),
        ("cost", 1.0, f64::INFINITY),
        ("cost", 1.0, -1.0),
    ] {
        assert!(matches!(
            Knob::new("bad-setting", 0, vec![ks("bad", error, cost)]),
            Err(PlanInputError::InvalidSettingValue {
                field: actual,
                ..
            }) if actual == field
        ));
    }

    let valid = Knob::new("valid", 0, vec![ks("base", 1.0, 1.0)]).expect("valid control knob");
    for (field, budget_s, error_target) in [
        ("budget_s", f64::NAN, 0.0),
        ("budget_s", f64::INFINITY, 0.0),
        ("budget_s", -1.0, 0.0),
        ("error_target", 1.0, f64::NAN),
        ("error_target", 1.0, f64::INFINITY),
        ("error_target", 1.0, -1.0),
    ] {
        assert!(matches!(
            allocate(&AllocProblem {
                knobs: vec![valid.clone()],
                budget_s,
                error_target,
            }),
            Err(AllocationError::InvalidInput(
                PlanInputError::InvalidProblemValue { field: actual, .. }
            )) if actual == field
        ));
    }

    let overflow = AllocProblem {
        knobs: vec![
            Knob::new("a", 0, vec![ks("base", f64::MAX, f64::MAX)])
                .expect("individual values are finite"),
            Knob::new("b", 0, vec![ks("base", f64::MAX, f64::MAX)])
                .expect("individual values are finite"),
        ],
        budget_s: f64::MAX,
        error_target: 0.0,
    };
    assert!(matches!(
        allocate(&overflow),
        Err(AllocationError::InvalidInput(
            PlanInputError::AggregateOverflow { .. }
        ))
    ));
    verdict(
        "pl-008-input-domains",
        true,
        "setting, budget, target, and aggregate domains reject explicitly",
    );
}

/// pl-009: track IDs are rejected before any storage or iteration can
/// be sized by the numeric ID. Validate both the constructor and the
/// defensive check for callers that construct the public struct.
#[test]
fn pl_009_huge_track_ids_refuse_before_sizing() {
    assert!(matches!(
        Knob::new(
            "too-many-tracks",
            MAX_EXECUTION_TRACKS,
            vec![ks("base", 1.0, 1.0)],
        ),
        Err(PlanInputError::TrackOutOfRange { .. })
    ));

    let forged = AllocProblem {
        knobs: vec![Knob {
            name: "forged".to_string(),
            track: usize::MAX,
            settings: vec![ks("base", 1.0, 1.0)],
        }],
        budget_s: 1.0,
        error_target: 1.0,
    };
    assert!(matches!(
        allocate(&forged),
        Err(AllocationError::InvalidInput(
            PlanInputError::TrackOutOfRange {
                track: usize::MAX,
                ..
            }
        ))
    ));
    assert!(matches!(
        oracle_min_error(&forged),
        Err(PlanInputError::TrackOutOfRange { .. })
    ));
    verdict(
        "pl-009-bounded-tracks",
        true,
        "usize::MAX track IDs refuse before track-indexed work",
    );
}

/// pl-010: invalid online evidence is a non-mutating refusal, and a
/// valid update that changes dominance leaves a valid Pareto ladder.
#[test]
fn pl_010_online_observations_are_transactional() {
    let knob = Knob::new(
        "mesh",
        0,
        vec![ks("coarse", 0.8, 1.0), ks("fine", 0.2, 4.0)],
    )
    .expect("valid mesh knob");
    let mut allocator = Allocator::new(AllocProblem {
        knobs: vec![knob],
        budget_s: 4.0,
        error_target: 0.2,
    })
    .expect("valid allocation problem");
    let before = allocator.problem().knobs[0].settings[0].error;
    assert!(matches!(
        allocator.observe_error(0, 0, f64::NAN),
        Err(PlanInputError::InvalidSettingValue { .. })
    ));
    assert_eq!(
        allocator.problem().knobs[0].settings[0].error.to_bits(),
        before.to_bits()
    );
    assert!(matches!(
        allocator.observe_error(1, 0, 0.1),
        Err(PlanInputError::KnobIndexOutOfRange { .. })
    ));

    allocator
        .observe_error(0, 0, 0.1)
        .expect("valid measurement commits");
    assert_eq!(allocator.problem().knobs[0].settings.len(), 1);
    let plan = allocator.replan().expect("re-pruned problem remains valid");
    assert_eq!(plan.choice, vec![0]);
    assert!(plan.wall_clock <= allocator.problem().budget_s);
    verdict(
        "pl-010-transactional-observation",
        true,
        "invalid evidence does not mutate; valid evidence re-prunes dominance",
    );
}

/// pl-011: no knobs is the intentional identity problem, not an
/// indexing accident: one empty choice has zero additive error and
/// zero tropical wall cost.
#[test]
fn pl_011_empty_problem_is_the_zero_cost_identity() {
    let problem = AllocProblem {
        knobs: Vec::new(),
        budget_s: 0.0,
        error_target: 1.0,
    };
    let plan = allocate(&problem).expect("identity problem is feasible");
    assert!(plan.choice.is_empty());
    assert_eq!(plan.total_error.to_bits(), 0.0f64.to_bits());
    assert_eq!(plan.wall_clock.to_bits(), 0.0f64.to_bits());
    assert_eq!(
        oracle_min_error(&problem)
            .expect("identity problem is valid")
            .expect("identity plan exists")
            .0,
        Vec::<usize>::new()
    );
    verdict(
        "pl-011-empty-identity",
        true,
        "empty knobs yield the unique zero-error, zero-wall plan",
    );
}

fn pareto_settings(count: usize) -> Vec<KnobSetting> {
    (0..count)
        .map(|index| {
            ks(
                &format!("s{index}"),
                (count - index) as f64,
                (index + 1) as f64,
            )
        })
        .collect()
}

/// pl-012: every planner-controlled work dimension is bounded, and
/// the exact fixture oracle refuses before Cartesian enumeration.
#[test]
fn pl_012_work_domains_are_bounded() {
    assert!(matches!(
        Knob::new("oversized", 0, pareto_settings(MAX_SETTINGS_PER_KNOB + 1),),
        Err(PlanInputError::TooManySettings { .. })
    ));

    let singleton = Knob::new("single", 0, vec![ks("base", 1.0, 1.0)]).expect("valid singleton");
    let too_many_knobs = AllocProblem {
        knobs: vec![singleton; MAX_ALLOCATION_KNOBS + 1],
        budget_s: 1.0,
        error_target: 1.0,
    };
    assert!(matches!(
        allocate(&too_many_knobs),
        Err(AllocationError::InvalidInput(
            PlanInputError::TooManyKnobs { .. }
        ))
    ));

    let dense = Knob::new("dense", 0, pareto_settings(MAX_SETTINGS_PER_KNOB))
        .expect("per-knob boundary is valid");
    let too_many_total = AllocProblem {
        knobs: vec![dense; MAX_TOTAL_SETTINGS / MAX_SETTINGS_PER_KNOB + 1],
        budget_s: f64::MAX,
        error_target: 0.0,
    };
    assert!(matches!(
        allocate(&too_many_total),
        Err(AllocationError::InvalidInput(
            PlanInputError::TooManyTotalSettings { .. }
        ))
    ));

    let oracle_problem = AllocProblem {
        knobs: (0..4)
            .map(|index| {
                Knob::new(&format!("oracle-{index}"), index, pareto_settings(33))
                    .expect("valid oracle knob")
            })
            .collect(),
        budget_s: 200.0,
        error_target: 0.0,
    };
    assert!(33usize.pow(4) > MAX_ORACLE_COMBINATIONS);
    assert!(matches!(
        oracle_min_error(&oracle_problem),
        Err(PlanInputError::OracleWorkLimitExceeded { .. })
    ));
    verdict(
        "pl-012-bounded-work",
        true,
        "knob, setting, aggregate, and exact-oracle work caps refuse explicitly",
    );
}

/// pl-013: the public evaluators validate their choice vectors rather
/// than indexing caller input directly.
#[test]
fn pl_013_public_evaluators_are_fallible() {
    let knobs =
        vec![Knob::new("mesh", 0, vec![ks("coarse", 1.0, 1.0)]).expect("valid evaluator knob")];
    assert!(matches!(
        fs_plan::alloc::plan_wall_clock(&knobs, &[]),
        Err(PlanInputError::ChoiceLengthMismatch { .. })
    ));
    assert!(matches!(
        fs_plan::alloc::plan_total_error(&knobs, &[1]),
        Err(PlanInputError::ChoiceIndexOutOfRange { .. })
    ));
    assert_eq!(
        fs_plan::alloc::plan_wall_clock(&knobs, &[0])
            .expect("valid choice")
            .to_bits(),
        1.0f64.to_bits()
    );
    verdict(
        "pl-013-fallible-evaluators",
        true,
        "choice length and indices refuse without panicking",
    );
}

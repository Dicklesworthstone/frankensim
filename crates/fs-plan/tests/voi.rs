//! VoI-query conformance (the knh1.6 bead; runs under `voi-queries`).
//! Acceptance: decision sensitivity from cached sweeps (near-free —
//! call counts prove it); rankings by flip-probability-per-dollar; the
//! probe menu unifies computational and physical experiments; the
//! ranking surfaces as the query hint and the probe scheduler; myopic
//! one-step only; the prospective-audit kill criterion demotes VoI
//! when recommendations stop outperforming.
#![cfg(feature = "voi-queries")]

use std::cell::Cell;

use fs_plan::voi::{
    AuditRecord, AuditVerdict, LiveDecision, MAX_VOI_EVALUATIONS, MAX_VOI_GRID, MAX_VOI_NAME_BYTES,
    MAX_VOI_NODES, MAX_VOI_PROBES, Probe, ProbeKind, RankedPurchase, UncertaintyNode, VoiError,
    audit_verdict, hint_for_query, rank_purchases, schedule_probes,
};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-plan/voi\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

/// The live decision: ship when margin = 2·drag_gap − mass_penalty +
/// 0.3·hazard > 0 (a cheap cached surrogate).
fn margin(v: &[f64]) -> f64 {
    2.0 * v[0] - v[1] + 0.3 * v[2]
}

fn nodes() -> Vec<UncertaintyNode> {
    vec![
        // Straddles the boundary: pivotal.
        UncertaintyNode {
            name: "drag-gap".to_string(),
            lo: -0.4,
            hi: 0.6,
            nominal: 0.1,
        },
        // Entirely on one side: irrelevant to the verdict.
        UncertaintyNode {
            name: "mass-penalty".to_string(),
            lo: 0.05,
            hi: 0.15,
            nominal: 0.1,
        },
        // Mildly pivotal.
        UncertaintyNode {
            name: "hazard".to_string(),
            lo: -1.0,
            hi: 1.0,
            nominal: 0.2,
        },
    ]
}

fn probe(name: &str, target: &str, cost: f64, shrink: f64) -> Probe {
    Probe {
        name: name.to_string(),
        target: target.to_string(),
        cost,
        shrink,
        kind: ProbeKind::Computational,
    }
}

fn ranked_probe(name: &str, cost: f64) -> RankedPurchase {
    let decision = LiveDecision {
        margin: &margin,
        arity: 3,
    };
    rank_purchases(
        &decision,
        &nodes(),
        &[probe(name, "drag-gap", cost, 0.01)],
        64,
    )
    .expect("valid ranked-probe fixture")
    .pop()
    .expect("one fixture purchase")
}

#[test]
fn voi_001_sensitivity_from_cached_sweeps() {
    let calls = Cell::new(0usize);
    let counting = |v: &[f64]| -> f64 {
        calls.set(calls.get() + 1);
        margin(v)
    };
    let decision = LiveDecision {
        margin: &counting,
        arity: 3,
    };
    let ns = nodes();
    let grid = 32;
    let p_pivotal = decision
        .flip_probability(&ns, 0, grid)
        .expect("valid pivotal sweep");
    let p_irrelevant = decision
        .flip_probability(&ns, 1, grid)
        .expect("valid irrelevant sweep");
    let p_mild = decision
        .flip_probability(&ns, 2, grid)
        .expect("valid mild sweep");
    println!(
        "{{\"metric\":\"sensitivity\",\"drag_gap\":{p_pivotal:.3},\
         \"mass_penalty\":{p_irrelevant:.3},\"hazard\":{p_mild:.3},\"surrogate_calls\":{}}}",
        calls.get()
    );
    assert!(
        p_pivotal > 0.3,
        "the straddling node is pivotal: {p_pivotal}"
    );
    assert_eq!(
        p_irrelevant.to_bits(),
        0.0f64.to_bits(),
        "the one-sided node cannot flip"
    );
    assert!(
        p_mild > 0.0 && p_mild < p_pivotal,
        "mildly pivotal: {p_mild}"
    );
    // NEAR-FREE: three sweeps cost 3·(grid + 1) surrogate calls — no
    // solver in the loop, no combinatorial sweep.
    assert!(
        calls.get() <= 3 * (grid + 1),
        "cached-sweep budget held: {} calls",
        calls.get()
    );
    verdict(
        "voi-001",
        "flip probabilities from cached interval sweeps: pivotal 0.4-class, one-sided \
         exactly 0, mild in between — at <= 3(grid+1) surrogate calls (near-free)",
    );
}

#[test]
fn voi_002_ranking_is_flip_prob_per_dollar() {
    let decision = LiveDecision {
        margin: &margin,
        arity: 3,
    };
    let ns = nodes();
    let menu = vec![
        // Cheap AND decisive: must rank first.
        Probe {
            name: "climb-rung-drag".to_string(),
            target: "drag-gap".to_string(),
            cost: 10.0,
            shrink: 0.2,
            kind: ProbeKind::Computational,
        },
        // Decisive but expensive: positive score, ranks below.
        Probe {
            name: "wind-tunnel-drag".to_string(),
            target: "drag-gap".to_string(),
            cost: 500.0,
            shrink: 0.05,
            kind: ProbeKind::Physical,
        },
        // Cheap but irrelevant: zero score, ranks last.
        Probe {
            name: "refine-mass-model".to_string(),
            target: "mass-penalty".to_string(),
            cost: 5.0,
            shrink: 0.1,
            kind: ProbeKind::Computational,
        },
    ];
    let ranked = rank_purchases(&decision, &ns, &menu, 64).expect("valid probe menu");
    let names: Vec<&str> = ranked.iter().map(|r| r.probe().name.as_str()).collect();
    println!(
        "{{\"metric\":\"ranking\",\"order\":{names:?},\"scores\":{:?}}}",
        ranked
            .iter()
            .map(|r| (r.score() * 1e4).round() / 1e4)
            .collect::<Vec<_>>()
    );
    assert_eq!(names[0], "climb-rung-drag", "cheap+decisive wins");
    assert_eq!(names[1], "wind-tunnel-drag", "decisive-but-pricey second");
    assert_eq!(names[2], "refine-mass-model", "irrelevant last");
    assert!(ranked[0].score() > ranked[1].score() && ranked[1].score() > 0.0);
    assert_eq!(
        ranked[2].score().to_bits(),
        0.0f64.to_bits(),
        "an irrelevant probe buys nothing"
    );
    verdict(
        "voi-002",
        "the ranking is flip-probability-per-dollar: cheap+decisive > decisive+pricey > \
         irrelevant (score exactly 0)",
    );
}

#[test]
fn voi_003_menu_unifies_compute_and_physical() {
    // A physical anchor that shrinks the pivotal node far harder can
    // BEAT a cheap but weak computational probe — the menu prices
    // evidence, not its substrate.
    let decision = LiveDecision {
        margin: &margin,
        arity: 3,
    };
    let ns = nodes();
    let menu = vec![
        Probe {
            name: "weak-rung-climb".to_string(),
            target: "drag-gap".to_string(),
            cost: 8.0,
            shrink: 0.95, // barely tightens
            kind: ProbeKind::Computational,
        },
        Probe {
            name: "wind-tunnel-anchor".to_string(),
            target: "drag-gap".to_string(),
            cost: 60.0,
            shrink: 0.05, // near-eliminates the interval
            kind: ProbeKind::Physical,
        },
    ];
    let ranked = rank_purchases(&decision, &ns, &menu, 64).expect("valid mixed menu");
    assert_eq!(ranked[0].probe().kind, ProbeKind::Physical);
    assert!(
        ranked[0].score() > ranked[1].score(),
        "the physical anchor wins on flip-prob-per-dollar: {:.4} vs {:.4}",
        ranked[0].score(),
        ranked[1].score()
    );
    verdict(
        "voi-003",
        "a decisive wind-tunnel anchor outranks a weak rung climb despite 7.5x the \
         price — computational and physical evidence priced on one menu",
    );
}

#[test]
fn voi_004_surfacing_hint_and_scheduler() {
    let decision = LiveDecision {
        margin: &margin,
        arity: 3,
    };
    let ns = nodes();
    let menu = vec![
        Probe {
            name: "climb-rung-drag".to_string(),
            target: "drag-gap".to_string(),
            cost: 10.0,
            shrink: 0.2,
            kind: ProbeKind::Computational,
        },
        Probe {
            name: "hazard-samples".to_string(),
            target: "hazard".to_string(),
            cost: 25.0,
            shrink: 0.3,
            kind: ProbeKind::Computational,
        },
        Probe {
            name: "wind-tunnel-drag".to_string(),
            target: "drag-gap".to_string(),
            cost: 500.0,
            shrink: 0.05,
            kind: ProbeKind::Physical,
        },
    ];
    let ranked = rank_purchases(&decision, &ns, &menu, 64).expect("valid scheduling menu");
    // (i) The query-result hint (the Proposal-8 anytime shape, now
    // decision-priced).
    let hint = hint_for_query(&ranked).expect("valid canonical hint");
    println!("{{\"metric\":\"hint\",\"text\":\"{hint}\"}}");
    assert!(hint.contains("climb-rung-drag") && hint.contains("$10"));
    assert!(hint.contains("flip-probability"));
    // (ii) The probe scheduler under a budget: greedy top-k affordable.
    let scheduled = schedule_probes(&ranked, 40.0).expect("valid finite schedule");
    let names: Vec<&str> = scheduled.iter().map(|r| r.probe().name.as_str()).collect();
    assert_eq!(
        names,
        vec!["climb-rung-drag", "hazard-samples"],
        "greedy affordable top-k under $40"
    );
    // Input order carries no authority: both surfaces reapply the
    // canonical score/cost/name comparator.
    let mut reordered = ranked.clone();
    reordered.reverse();
    assert_eq!(hint_for_query(&reordered).expect("reordered hint"), hint);
    assert_eq!(
        schedule_probes(&reordered, 40.0).expect("reordered schedule"),
        scheduled
    );
    // Myopic-only is structural: rank_purchases takes ONE state and
    // returns ONE ranked step — there is no sequential-tree API to
    // misuse.
    verdict(
        "voi-004",
        "the top purchase surfaces as the priced query hint; the scheduler buys the \
         greedy affordable top-k under budget; the API is one-step myopic by shape",
    );
}

#[test]
fn voi_005_prospective_audit_kill_criterion() {
    // Recommended purchases realize more decision changes: keep.
    let winning: Vec<AuditRecord> = (0..10)
        .map(|k| AuditRecord {
            recommended_changed_decision: k < 6,
            alternative_changed_decision: k < 3,
        })
        .collect();
    assert_eq!(audit_verdict(&winning), AuditVerdict::KeepScheduling);
    // They stop outperforming at matched cost: DEMOTE to reporting.
    let losing: Vec<AuditRecord> = (0..10)
        .map(|k| AuditRecord {
            recommended_changed_decision: k < 3,
            alternative_changed_decision: k < 4,
        })
        .collect();
    assert_eq!(audit_verdict(&losing), AuditVerdict::DemoteToReporting);
    // No evidence, no authority.
    assert_eq!(audit_verdict(&[]), AuditVerdict::DemoteToReporting);
    verdict(
        "voi-005",
        "the prospective audit keeps VoI's scheduling authority only while recommended \
         purchases measurably outperform agent-chosen alternatives — and with no audit \
         evidence there is no authority",
    );
}

#[test]
fn voi_006_node_name_and_grid_boundaries() {
    let boundary_nodes: Vec<UncertaintyNode> = (0..MAX_VOI_NODES)
        .map(|index| UncertaintyNode {
            name: format!("node-{index}"),
            lo: 0.0,
            hi: 1.0,
            nominal: 0.5,
        })
        .collect();
    let wide_margin = |values: &[f64]| values[0] - 0.5;
    let boundary_decision = LiveDecision {
        margin: &wide_margin,
        arity: MAX_VOI_NODES,
    };
    assert!(
        boundary_decision.nominal_verdict(&boundary_nodes).is_ok(),
        "the exact node-count boundary is admitted"
    );
    let mut too_many = boundary_nodes.clone();
    too_many.push(UncertaintyNode {
        name: "one-too-many".to_string(),
        lo: 0.0,
        hi: 1.0,
        nominal: 0.5,
    });
    let oversized_decision = LiveDecision {
        margin: &wide_margin,
        arity: too_many.len(),
    };
    assert!(matches!(
        oversized_decision.nominal_verdict(&too_many),
        Err(VoiError::SizeLimit { .. })
    ));

    let exact_name = vec![UncertaintyNode {
        name: "n".repeat(MAX_VOI_NAME_BYTES),
        lo: 0.0,
        hi: 1.0,
        nominal: 0.5,
    }];
    let one = LiveDecision {
        margin: &wide_margin,
        arity: 1,
    };
    assert!(one.nominal_verdict(&exact_name).is_ok());
    let mut long_name = exact_name;
    long_name[0].name.push('n');
    assert!(matches!(
        one.nominal_verdict(&long_name),
        Err(VoiError::InvalidName { .. })
    ));

    let single = vec![UncertaintyNode {
        name: "x".to_string(),
        lo: -1.0,
        hi: 1.0,
        nominal: 0.0,
    }];
    assert!(
        one.flip_probability(&single, 0, MAX_VOI_GRID).is_ok(),
        "the exact grid boundary is admitted"
    );
    for grid in [0, MAX_VOI_GRID + 1] {
        assert!(matches!(
            one.flip_probability(&single, 0, grid),
            Err(VoiError::InvalidGrid { .. })
        ));
    }
    verdict(
        "voi-006",
        "node count, UTF-8 name bytes, and sweep grid admit exact boundaries and refuse limit+1",
    );
}

#[test]
fn voi_007_menu_and_aggregate_evaluation_boundaries() {
    let calls = Cell::new(0usize);
    let counting = |values: &[f64]| {
        calls.set(calls.get() + 1);
        values[0]
    };
    let decision = LiveDecision {
        margin: &counting,
        arity: 1,
    };
    let ns = vec![UncertaintyNode {
        name: "x".to_string(),
        lo: -1.0,
        hi: 1.0,
        nominal: 0.0,
    }];
    let menu: Vec<Probe> = (0..MAX_VOI_PROBES)
        .map(|index| probe(&format!("probe-{index}"), "x", index as f64 + 1.0, 0.5))
        .collect();
    assert_eq!(MAX_VOI_PROBES * 2 * (1 + 1), MAX_VOI_EVALUATIONS);
    let ranked = rank_purchases(&decision, &ns, &menu, 1)
        .expect("exact menu/evaluation boundary is admitted");
    assert_eq!(ranked.len(), MAX_VOI_PROBES);
    assert_eq!(calls.get(), MAX_VOI_EVALUATIONS);

    calls.set(0);
    assert!(matches!(
        rank_purchases(&decision, &ns, &menu, 2),
        Err(VoiError::EvaluationLimitExceeded { .. })
    ));
    assert_eq!(calls.get(), 0, "aggregate refusal precedes callbacks");

    let mut oversized = menu;
    oversized.push(probe("one-too-many", "x", 1.0, 0.5));
    assert!(matches!(
        rank_purchases(&decision, &ns, &oversized, 1),
        Err(VoiError::SizeLimit { .. })
    ));
    assert_eq!(calls.get(), 0, "menu-size refusal precedes callbacks");
    verdict(
        "voi-007",
        "menu and aggregate evaluation caps admit exact limits and refuse before callback at limit+1",
    );
}

#[test]
fn voi_008_decision_interval_and_callback_refusals() {
    let calls = Cell::new(0usize);
    let counting = |values: &[f64]| {
        calls.set(calls.get() + 1);
        values.first().copied().unwrap_or(0.0)
    };
    let decision = LiveDecision {
        margin: &counting,
        arity: 1,
    };
    assert!(matches!(
        decision.nominal_verdict(&[]),
        Err(VoiError::SizeLimit { .. })
    ));
    assert_eq!(calls.get(), 0);
    let ns = vec![UncertaintyNode {
        name: "x".to_string(),
        lo: 0.0,
        hi: 1.0,
        nominal: 0.5,
    }];
    let wrong_arity = LiveDecision {
        margin: &counting,
        arity: 2,
    };
    assert!(matches!(
        wrong_arity.nominal_verdict(&ns),
        Err(VoiError::ArityMismatch { .. })
    ));
    assert!(matches!(
        decision.flip_probability(&ns, 1, 4),
        Err(VoiError::NodeIndexOutOfRange { .. })
    ));

    let mut duplicate = ns.clone();
    duplicate.push(ns[0].clone());
    let duplicate_decision = LiveDecision {
        margin: &counting,
        arity: 2,
    };
    assert!(matches!(
        duplicate_decision.nominal_verdict(&duplicate),
        Err(VoiError::DuplicateName { .. })
    ));
    for (lo, nominal, hi) in [
        (f64::NAN, 0.5, 1.0),
        (1.0, 0.5, 0.0),
        (0.0, 2.0, 1.0),
        (-f64::MAX, 0.0, f64::MAX),
    ] {
        let invalid = vec![UncertaintyNode {
            name: "x".to_string(),
            lo,
            hi,
            nominal,
        }];
        assert!(matches!(
            decision.nominal_verdict(&invalid),
            Err(VoiError::InvalidInterval { .. })
        ));
    }

    let nan_margin = |_: &[f64]| f64::NAN;
    let malformed = LiveDecision {
        margin: &nan_margin,
        arity: 1,
    };
    assert!(matches!(
        malformed.nominal_verdict(&ns),
        Err(VoiError::NonFiniteMargin { .. })
    ));
    let interior_nan = |values: &[f64]| {
        if values[0] > 0.5 { f64::INFINITY } else { 0.0 }
    };
    let malformed = LiveDecision {
        margin: &interior_nan,
        arity: 1,
    };
    assert!(matches!(
        malformed.flip_probability(&ns, 0, 2),
        Err(VoiError::NonFiniteMargin { .. })
    ));
    verdict(
        "voi-008",
        "arity, identity, intervals, indices, duplicate nodes, and nonfinite callbacks refuse structurally",
    );
}

#[test]
fn voi_009_menu_targets_and_probe_economics_refuse_before_sweep() {
    let calls = Cell::new(0usize);
    let counting = |values: &[f64]| {
        calls.set(calls.get() + 1);
        values[0]
    };
    let decision = LiveDecision {
        margin: &counting,
        arity: 1,
    };
    let ns = vec![UncertaintyNode {
        name: "x".to_string(),
        lo: -1.0,
        hi: 1.0,
        nominal: 0.0,
    }];
    assert!(matches!(
        rank_purchases(&decision, &ns, &[], 4),
        Err(VoiError::SizeLimit { .. })
    ));
    let duplicate_nodes = vec![ns[0].clone(), ns[0].clone()];
    let duplicate_target_decision = LiveDecision {
        margin: &counting,
        arity: 2,
    };
    assert!(matches!(
        rank_purchases(
            &duplicate_target_decision,
            &duplicate_nodes,
            &[probe("ambiguous", "x", 1.0, 0.5)],
            4,
        ),
        Err(VoiError::DuplicateName { .. })
    ));
    let unknown = vec![probe("unknown", "missing", 1.0, 0.5)];
    let unknown_before = unknown.clone();
    assert!(matches!(
        rank_purchases(&decision, &ns, &unknown, 4),
        Err(VoiError::TargetResolution { matches: 0, .. })
    ));
    assert_eq!(unknown, unknown_before, "refusal does not mutate the menu");
    assert_eq!(calls.get(), 0, "unknown target refuses before sweep");

    let duplicated = vec![probe("same", "x", 1.0, 0.5), probe("same", "x", 2.0, 0.5)];
    assert!(matches!(
        rank_purchases(&decision, &ns, &duplicated, 4),
        Err(VoiError::DuplicateName { .. })
    ));
    for cost in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            rank_purchases(&decision, &ns, &[probe("bad-cost", "x", cost, 0.5)], 4),
            Err(VoiError::InvalidProbeValue { field: "cost", .. })
        ));
    }
    for shrink in [0.0, 1.0, -0.1, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            rank_purchases(&decision, &ns, &[probe("bad-shrink", "x", 1.0, shrink)], 4),
            Err(VoiError::InvalidProbeValue {
                field: "shrink",
                ..
            })
        ));
    }
    for bad_name in ["", " padded", "trailing "] {
        assert!(matches!(
            rank_purchases(&decision, &ns, &[probe(bad_name, "x", 1.0, 0.5)], 4),
            Err(VoiError::InvalidName { .. })
        ));
    }
    assert_eq!(calls.get(), 0, "all malformed menus refuse before sweep");
    verdict(
        "voi-009",
        "unknown/ambiguous targets, duplicate identities, and invalid probe economics refuse before evaluation",
    );
}

#[test]
fn voi_010_probe_name_and_score_arithmetic_boundaries() {
    let calls = Cell::new(0usize);
    let counting = |values: &[f64]| {
        calls.set(calls.get() + 1);
        values[0]
    };
    let decision = LiveDecision {
        margin: &counting,
        arity: 1,
    };
    let ns = vec![UncertaintyNode {
        name: "x".to_string(),
        lo: -1.0,
        hi: 1.0,
        nominal: 0.0,
    }];
    let exact = "p".repeat(MAX_VOI_NAME_BYTES);
    assert!(
        rank_purchases(&decision, &ns, &[probe(&exact, "x", 1.0, 0.5)], 4).is_ok(),
        "the exact probe-name byte boundary is admitted"
    );
    let before_long_name = calls.get();
    let too_long = "p".repeat(MAX_VOI_NAME_BYTES + 1);
    assert!(matches!(
        rank_purchases(&decision, &ns, &[probe(&too_long, "x", 1.0, 0.5)], 4),
        Err(VoiError::InvalidName { .. })
    ));
    assert_eq!(calls.get(), before_long_name);
    let too_long_target = "t".repeat(MAX_VOI_NAME_BYTES + 1);
    assert!(matches!(
        rank_purchases(
            &decision,
            &ns,
            &[probe("long-target", &too_long_target, 1.0, 0.5)],
            4,
        ),
        Err(VoiError::InvalidName { .. })
    ));

    let decision = LiveDecision {
        margin: &margin,
        arity: 3,
    };
    assert!(matches!(
        rank_purchases(
            &decision,
            &nodes(),
            &[probe("underflow-cost", "drag-gap", f64::from_bits(1), 0.01,)],
            64,
        ),
        Err(VoiError::ArithmeticRefusal { .. })
    ));
    verdict(
        "voi-010",
        "probe names admit the exact byte bound, limit+1 refuses before callback, and derived score overflow refuses",
    );
}

#[test]
fn voi_011_scheduler_is_transactional_and_budget_monotone() {
    let valid = vec![ranked_probe("a", 10.0), ranked_probe("b", 5.0)];
    assert!(
        schedule_probes(&valid, 0.0)
            .expect("zero budget is valid")
            .is_empty()
    );
    for budget in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            schedule_probes(&valid, budget),
            Err(VoiError::InvalidBudget { .. })
        ));
    }
    assert!(matches!(
        schedule_probes(&[], 1.0),
        Err(VoiError::SizeLimit { .. })
    ));

    let duplicated = vec![ranked_probe("same", 10.0), ranked_probe("same", 10.0)];
    let before = duplicated.clone();
    assert!(matches!(
        schedule_probes(&duplicated, 100.0),
        Err(VoiError::DuplicateRankedProbe { .. })
    ));
    assert_eq!(
        duplicated, before,
        "duplicate refusal does not mutate input"
    );

    assert!(matches!(
        hint_for_query(&[]),
        Err(VoiError::SizeLimit { .. })
    ));
    assert!(matches!(
        hint_for_query(&duplicated),
        Err(VoiError::DuplicateRankedProbe { .. })
    ));
    let no_progress = vec![ranked_probe("tiny", 1.0)];
    assert!(matches!(
        schedule_probes(&no_progress, f64::MAX),
        Err(VoiError::ArithmeticRefusal { .. })
    ));

    let scheduled = schedule_probes(&valid, 15.0).expect("exact finite budget");
    assert_eq!(scheduled.len(), 2);
    assert_eq!(
        scheduled
            .iter()
            .map(|purchase| purchase.probe().name.as_str())
            .collect::<Vec<_>>(),
        vec!["b", "a"]
    );
    verdict(
        "voi-011",
        "invalid budgets and duplicate identities refuse transactionally; hints and schedules canonicalize authority",
    );
}

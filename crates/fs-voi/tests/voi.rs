//! Battery for value-of-information (fs-voi). Covers ranking-flip probability,
//! the FULL multi-alternative expected opportunity loss (robust → 0, close →
//! positive, high-variance third alternative → NOT robust), the top-two
//! surrogate's honest demotion (bead sj31i.5), the decision-vs-estimator
//! distinction (info on a decision-irrelevant design is worthless),
//! STOP-when-robust, VOI beating the uncertainty-proportional baseline,
//! model-form escalation, overflow-safe uncertainty composition, and a
//! deterministic Monte Carlo oracle comparison.

use fs_voi::{
    Action, ActionKind, Component, DesignEstimate, EVPI_SEMANTICS_VERSION, Recommendation,
    Uncertainty, action_value, decision_posture, expected_opportunity_loss, heuristic_choice,
    ranking_flip_probability, recommend, top_two_evpi_surrogate,
};

fn unc(n: f64, s: f64, m: f64) -> Uncertainty {
    Uncertainty {
        numerical: n,
        statistical: s,
        model: m,
    }
}
fn design(name: &str, mean: f64, u: Uncertainty) -> DesignEstimate {
    DesignEstimate::new(name, mean, u)
}
fn act(name: &str, kind: ActionKind, target: &str, cost: f64) -> Action {
    Action {
        name: name.into(),
        kind,
        target_design: target.into(),
        reduction: 0.9,
        cost,
    }
}

#[test]
fn evpi_semantics_version_is_locked() {
    // v2 (bead sj31i.5): full multi-alternative expected opportunity
    // loss carries robustness; the top-two closed form is a renamed
    // surrogate; uncertainty composition is a scaled norm.
    assert_eq!(EVPI_SEMANTICS_VERSION, 2);
}

#[test]
fn ranking_flip_probability_reflects_separation() {
    let a = design("a", 0.0, unc(1.0, 0.0, 0.0));
    let far = design("far", 3.0, unc(1.0, 0.0, 0.0));
    // 3σ-ish apart -> flip is unlikely.
    assert!(ranking_flip_probability(&a, &far) < 0.05);
    // a near tie -> flip probability near 0.5.
    let near = design("near", 0.1, unc(1.0, 0.0, 0.0));
    let p = ranking_flip_probability(&a, &near);
    assert!(p > 0.4 && p < 0.5);
}

#[test]
fn opportunity_loss_is_zero_for_a_robust_decision_and_positive_when_close() {
    let robust = [
        design("a", 0.0, unc(0.05, 0.0, 0.0)),
        design("b", 20.0, unc(0.05, 0.0, 0.0)),
    ];
    assert!(expected_opportunity_loss(&robust) < 1e-6);
    let close = [
        design("a", 0.0, unc(0.1, 0.1, 0.1)),
        design("c", 0.15, unc(0.1, 0.1, 0.1)),
    ];
    assert!(expected_opportunity_loss(&close) > 0.0);
    // the posture names the two closest designs.
    assert_eq!(
        decision_posture(&close).map(|p| (p.best, p.runner_up)),
        Some(("a".to_string(), "c".to_string()))
    );
}

/// The bead's core falsifier: a third alternative with a WORSE mean but
/// a LARGE uncertainty holds material probability of being optimal. The
/// top-two surrogate reads ~zero (mean ordering is not stochastic
/// dominance); the full evaluator must expose the loss and the
/// recommender must refuse to declare robustness.
#[test]
fn high_variance_third_alternative_defeats_the_top_two_surrogate() {
    let designs = [
        design("a", 0.0, unc(0.001, 0.0, 0.0)),
        design("b", 0.08, unc(0.001, 0.0, 0.0)),
        design("wild", 3.0, unc(5.0, 0.0, 0.0)),
    ];
    let surrogate = top_two_evpi_surrogate(&designs);
    let full = expected_opportunity_loss(&designs);
    assert!(
        surrogate < 1e-6,
        "top-two pair is 80σ apart; surrogate should read ~0, got {surrogate}"
    );
    // Closed form for the two-design (a, wild) reduction:
    // E[(0 − X_wild)⁺] = σφ(δ/σ) + δΦ(δ/σ) with δ = −3, σ ≈ 5.
    assert!(
        full > 0.5,
        "the wild alternative is optimal with P≈0.27; full EOL must be material, got {full}"
    );
    let rec = recommend(
        &designs,
        &[act("test-wild", ActionKind::Test, "wild", 1.0)],
        0.01,
    );
    assert!(
        !matches!(&rec, Recommendation::Stop { reason } if reason.contains("robust")),
        "global robustness must be impossible while an included alternative \
         holds material optimality probability: {rec:?}"
    );
}

/// With exactly two designs the full quadrature must agree with the
/// closed-form pairwise opportunity loss to quadrature resolution.
#[test]
fn full_evaluator_matches_the_closed_form_on_two_designs() {
    for (mean_b, std_a, std_b) in [
        (0.15, 0.17, 0.17),
        (0.5, 1.0, 2.0),
        (3.0, 0.25, 4.0),
        (0.0, 1.0, 1.0), // exact tie
    ] {
        let designs = [
            design("a", 0.0, unc(std_a, 0.0, 0.0)),
            design("b", mean_b, unc(std_b, 0.0, 0.0)),
        ];
        let full = expected_opportunity_loss(&designs);
        let closed = top_two_evpi_surrogate(&designs);
        let tolerance = 1e-6 * closed.abs().max(1.0);
        assert!(
            (full - closed).abs() <= tolerance,
            "two-design full ({full}) vs closed form ({closed}) at \
             (mean_b={mean_b}, std_a={std_a}, std_b={std_b})"
        );
    }
}

/// More alternatives can only add opportunity loss: the full evaluator
/// dominates the top-two surrogate on every menu.
#[test]
fn full_opportunity_loss_dominates_the_surrogate() {
    let menus: [&[DesignEstimate]; 3] = [
        &[
            design("a", 0.0, unc(0.1, 0.1, 0.1)),
            design("b", 0.15, unc(0.1, 0.1, 0.1)),
            design("c", 0.3, unc(0.5, 0.0, 0.0)),
        ],
        &[
            design("a", 1.0, unc(0.3, 0.0, 0.0)),
            design("b", 1.1, unc(0.4, 0.0, 0.0)),
            design("c", 1.2, unc(0.5, 0.0, 0.0)),
            design("d", 1.3, unc(0.6, 0.0, 0.0)),
        ],
        &[
            design("t1", 5.0, unc(1.0, 0.0, 0.0)),
            design("t2", 5.0, unc(1.0, 0.0, 0.0)),
            design("t3", 5.0, unc(1.0, 0.0, 0.0)),
        ],
    ];
    for menu in menus {
        let full = expected_opportunity_loss(menu);
        let surrogate = top_two_evpi_surrogate(menu);
        assert!(
            full >= surrogate - 1e-9,
            "full EOL ({full}) must dominate the top-two surrogate ({surrogate})"
        );
    }
    // Three exact ties: E[max of 3 std normals] ≈ 0.8463σ, strictly
    // above the two-design tie value √2·φ(0) ≈ 0.5642σ.
    let ties = &menus[2];
    let full = expected_opportunity_loss(ties);
    assert!((full - 0.8463).abs() < 5e-3, "three-way tie EOL: {full}");
}

/// Menu permutation cannot change the decision (the evaluator includes
/// everyone); values agree to floating-point re-association tolerance.
/// Canonical-order callers (the SensorForge menu) get bitwise identity
/// from their own ordering.
#[test]
fn full_opportunity_loss_is_permutation_invariant() {
    let base = [
        design("a", 0.0, unc(0.1, 0.0, 0.0)),
        design("b", 0.2, unc(0.3, 0.0, 0.0)),
        design("c", 1.0, unc(2.0, 0.0, 0.0)),
        design("d", 4.0, unc(6.0, 0.0, 0.0)),
    ];
    let reference = expected_opportunity_loss(&base);
    let mut permuted = base.clone();
    permuted.reverse();
    permuted.swap(1, 2);
    let value = expected_opportunity_loss(&permuted);
    assert!(
        (value - reference).abs() <= 1e-12 * reference.abs().max(1.0),
        "permutation moved full EOL: {reference} vs {value}"
    );
}

/// A soundly dominated alternative (far mean, ordinary uncertainty)
/// contributes nothing: including it is free, which is exactly why
/// nothing needs to be excluded.
#[test]
fn dominated_alternatives_contribute_nothing() {
    let core = [
        design("a", 0.0, unc(0.1, 0.1, 0.1)),
        design("b", 0.15, unc(0.1, 0.1, 0.1)),
    ];
    let with_dominated = [
        core[0].clone(),
        core[1].clone(),
        design("dominated", 1.0e6, unc(1.0, 0.0, 0.0)),
    ];
    let without = expected_opportunity_loss(&core);
    let with = expected_opportunity_loss(&with_dominated);
    assert!(
        (with - without).abs() <= 1e-9 * without.max(1.0),
        "a 1e6σ-dominated design changed the loss: {without} vs {with}"
    );
}

/// Deterministic Monte Carlo oracle: a seeded LCG + Box-Muller estimate
/// of E[X_best − min_j X_j] must agree with the quadrature within its
/// own confidence bound.
#[test]
fn monte_carlo_oracle_confirms_the_quadrature() {
    let designs = [
        design("a", 0.0, unc(0.5, 0.0, 0.0)),
        design("b", 0.3, unc(0.8, 0.0, 0.0)),
        design("c", 1.0, unc(2.5, 0.0, 0.0)),
        design("d", 2.0, unc(0.2, 0.0, 0.0)),
    ];
    let quadrature = expected_opportunity_loss(&designs);

    let mut state = 0x5d2f_e305_ce90_06fbu64; // deterministic seed
    let mut next_uniform = move || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        // Top 53 bits → (0, 1]; never 0 so ln() is finite.
        (((state >> 11) as f64) + 1.0) / (9007199254740992.0 + 1.0)
    };
    let samples = 200_000usize;
    let mut sum = 0.0f64;
    let mut sum_squares = 0.0f64;
    for _ in 0..samples {
        // Box-Muller pairs; draw enough normals for all four designs.
        let mut normals = [0.0f64; 4];
        for pair in 0..2 {
            let (u1, u2) = (next_uniform(), next_uniform());
            let radius = (-2.0 * u1.ln()).sqrt();
            let angle = std::f64::consts::TAU * u2;
            normals[2 * pair] = radius * angle.cos();
            normals[2 * pair + 1] = radius * angle.sin();
        }
        // The opportunity loss of committing to the apparent best
        // (design "a", the lowest mean): its SAMPLED value minus the
        // sampled menu minimum — nonnegative pointwise because the
        // minimum includes the chosen design, with expectation
        // `μ_a − E[min]`, exactly the evaluator's claim.
        let chosen = normals[0].mul_add(designs[0].uncertainty.total_std(), designs[0].mean);
        let minimum = designs
            .iter()
            .zip(normals)
            .map(|(d, z)| z.mul_add(d.uncertainty.total_std(), d.mean))
            .fold(f64::INFINITY, f64::min);
        let loss = chosen - minimum;
        sum += loss;
        sum_squares += loss * loss;
    }
    let n = samples as f64;
    let mc_mean = sum / n;
    let variance = (sum_squares / n - mc_mean * mc_mean).max(0.0);
    let standard_error = (variance / n).sqrt();
    assert!(
        (quadrature - mc_mean).abs() <= 5.0 * standard_error + 1e-4,
        "quadrature {quadrature} vs Monte Carlo {mc_mean} ± {standard_error}"
    );
}

/// The bead's numeric falsifier: finite components near √MAX must
/// compose to their representable total instead of overflowing the
/// naive variance sum, end to end through the evaluators.
#[test]
fn near_sqrt_max_uncertainties_stay_finite_end_to_end() {
    let huge = 1.0e154; // (1e154)² = 1e308 — the naive square sum overflows
    let composed = unc(huge, huge, huge).total_std();
    assert!(
        composed.is_finite(),
        "scaled norm must survive near-√MAX components"
    );
    assert!((composed / huge - 3.0f64.sqrt()).abs() < 1e-12);

    let designs = [
        design("a", 0.0, unc(huge, 0.0, 0.0)),
        design("b", 1.0e153, unc(huge, huge, 0.0)),
    ];
    let surrogate = top_two_evpi_surrogate(&designs);
    assert!(
        surrogate.is_finite() && surrogate > 0.0,
        "pairwise deviation composition must not overflow: {surrogate}"
    );
    let full = expected_opportunity_loss(&designs);
    assert!(
        full.is_finite() && full > 0.0,
        "full evaluator must survive astronomically scaled menus: {full}"
    );
    let tolerance = 1e-6 * surrogate;
    assert!(
        (full - surrogate).abs() <= tolerance,
        "two-design agreement must hold at scale: {full} vs {surrogate}"
    );
}

/// Power-of-two rescaling is exact in every quadrature intermediate, so
/// the full evaluator is bitwise scale-equivariant — a strong
/// metamorphic determinism check across ~120 orders of magnitude.
#[test]
fn full_opportunity_loss_is_bitwise_power_of_two_scale_equivariant() {
    let scale = 2.0f64.powi(400);
    let base = [
        design("a", 0.0, unc(0.5, 0.0, 0.0)),
        design("b", 0.3, unc(0.8, 0.0, 0.0)),
        design("c", 1.0, unc(2.5, 0.0, 0.0)),
    ];
    let scaled: Vec<DesignEstimate> = base
        .iter()
        .map(|d| {
            design(
                &d.name,
                d.mean * scale,
                unc(
                    d.uncertainty.numerical * scale,
                    d.uncertainty.statistical * scale,
                    d.uncertainty.model * scale,
                ),
            )
        })
        .collect();
    let reference = expected_opportunity_loss(&base);
    let rescaled = expected_opportunity_loss(&scaled) / scale;
    assert_eq!(
        reference.to_bits(),
        rescaled.to_bits(),
        "power-of-two scaling must be exact: {reference} vs {rescaled}"
    );
}

/// Subnormal-scale menus stay finite and non-negative.
#[test]
fn subnormal_scale_menus_are_admissible() {
    let tiny = 1.0e-310; // subnormal
    let designs = [
        design("a", 0.0, unc(tiny, 0.0, 0.0)),
        design("b", 2.0e-310, unc(tiny, 0.0, 0.0)),
    ];
    let full = expected_opportunity_loss(&designs);
    assert!(full.is_finite() && full >= 0.0, "subnormal EOL: {full}");
}

#[test]
fn non_finite_means_do_not_win_the_decision_boundary() {
    let designs = [
        design("nan", f64::NAN, unc(100.0, 0.0, 0.0)),
        design("a", 0.0, unc(0.1, 0.1, 0.1)),
        design("c", 0.15, unc(0.1, 0.1, 0.1)),
    ];
    assert_eq!(
        decision_posture(&designs).map(|p| (p.best, p.runner_up)),
        Some(("a".to_string(), "c".to_string()))
    );
    assert!(expected_opportunity_loss(&designs).is_finite());
    assert!(top_two_evpi_surrogate(&designs).is_finite());
    let insufficient = [
        design("nan", f64::NAN, unc(100.0, 0.0, 0.0)),
        design("a", 0.0, unc(0.1, 0.1, 0.1)),
    ];
    assert_eq!(decision_posture(&insufficient), None);
    assert!(expected_opportunity_loss(&insufficient).abs() <= f64::EPSILON);
    assert!(top_two_evpi_surrogate(&insufficient).abs() <= f64::EPSILON);
}

#[test]
fn information_on_a_decision_irrelevant_design_is_worthless() {
    let designs = [
        design("a", 0.0, unc(0.1, 0.1, 0.1)),  // best, in the decision
        design("c", 0.15, unc(0.1, 0.1, 0.1)), // close runner-up, in the decision
        design("b", 20.0, unc(5.0, 0.0, 0.0)), // clearly last but VERY uncertain
    ];
    // reducing the boundary design's uncertainty has value...
    let on_boundary = action_value(&designs, &act("refine-a", ActionKind::Refine, "a", 1.0));
    assert!(on_boundary.value > 0.0);
    // ...but reducing the uncertain-yet-irrelevant design's has ~none.
    let on_irrelevant = action_value(&designs, &act("refine-b", ActionKind::Refine, "b", 1.0));
    assert!(on_irrelevant.value < 1e-9);
}

#[test]
fn a_robust_decision_recommends_stop() {
    let robust = [
        design("a", 0.0, unc(0.05, 0.0, 0.0)),
        design("b", 20.0, unc(0.05, 0.0, 0.0)),
    ];
    let actions = [act("refine-a", ActionKind::Refine, "a", 1.0)];
    assert!(matches!(
        recommend(&robust, &actions, 1e-3),
        Recommendation::Stop { .. }
    ));
}

#[test]
fn voi_beats_the_uncertainty_proportional_baseline() {
    let designs = [
        design("a", 0.0, unc(0.1, 0.1, 0.1)),
        design("c", 0.15, unc(0.1, 0.1, 0.1)),
        design("b", 20.0, unc(5.0, 0.0, 0.0)), // most uncertain, decision-irrelevant
    ];
    let actions = [
        act("refine-a", ActionKind::Refine, "a", 1.0),
        act("refine-c", ActionKind::Refine, "c", 1.0),
        act("refine-b", ActionKind::Refine, "b", 1.0),
    ];
    // the uncertainty-proportional baseline chases the most-uncertain design (b).
    assert_eq!(
        heuristic_choice(&designs, &actions).unwrap().name,
        "refine-b"
    );
    // VOI spends on the DECISION boundary (a or c) instead.
    let rec = recommend(&designs, &actions, 1e-4);
    assert!(
        matches!(&rec, Recommendation::Act { action, .. } if action == "refine-a" || action == "refine-c"),
        "{rec:?}"
    );
}

#[test]
fn zero_cost_decision_changing_action_wins_per_cost_ranking() {
    let designs = [
        design("a", 0.0, unc(0.1, 0.1, 0.1)),
        design("c", 0.15, unc(0.1, 0.1, 0.1)),
    ];
    let free_action = act("free-refine-c", ActionKind::Refine, "c", 0.0);
    let free = action_value(&designs, &free_action);
    let actions = [
        act("paid-refine-a", ActionKind::Refine, "a", 1.0),
        free_action,
    ];
    assert!(free.value > 0.0);
    assert!(free.value_per_cost.is_infinite());
    assert!(matches!(
        recommend(&designs, &actions, 1e-4),
        Recommendation::Act { action, .. } if action == "free-refine-c"
    ));
}

#[test]
fn invalid_cost_actions_are_not_recommended() {
    let designs = [
        design("a", 0.0, unc(0.1, 0.1, 0.1)),
        design("c", 0.15, unc(0.1, 0.1, 0.1)),
    ];
    let negative = act("negative-cost", ActionKind::Refine, "c", -1.0);
    let nan_cost = act("nan-cost", ActionKind::Refine, "c", f64::NAN);
    assert!(action_value(&designs, &negative).value > 0.0);
    assert!(action_value(&designs, &negative).value_per_cost <= f64::EPSILON);
    assert!(action_value(&designs, &nan_cost).value > 0.0);
    assert!(action_value(&designs, &nan_cost).value_per_cost <= f64::EPSILON);
    assert!(matches!(
        recommend(&designs, &[negative, nan_cost], 1e-4),
        Recommendation::Stop { .. }
    ));
}

#[test]
fn voi_escalates_model_fidelity_when_model_uncertainty_dominates() {
    // the decision boundary is blocked by MODEL uncertainty (0.5), not statistical.
    let designs = [
        design("a", 0.0, unc(0.01, 0.01, 0.5)),
        design("c", 0.1, unc(0.01, 0.01, 0.5)),
    ];
    assert_eq!(designs[1].uncertainty.dominant(), Component::Model);
    let actions = [
        act("sample-c", ActionKind::Sample, "c", 1.0), // reduces statistical (tiny)
        act("test-c", ActionKind::Test, "c", 1.0),     // reduces model (dominant)
    ];
    // sampling barely helps; a physical test resolves the decision -> VOI escalates to Test.
    assert!(action_value(&designs, &actions[1]).value > action_value(&designs, &actions[0]).value);
    let rec = recommend(&designs, &actions, 1e-4);
    assert!(
        matches!(&rec, Recommendation::Act { action, .. } if action == "test-c"),
        "{rec:?}"
    );
}

#[test]
fn voi_is_deterministic() {
    let designs = [
        design("a", 0.0, unc(0.1, 0.1, 0.1)),
        design("c", 0.12, unc(0.1, 0.1, 0.1)),
    ];
    let actions = [act("refine-a", ActionKind::Refine, "a", 1.0)];
    assert_eq!(
        recommend(&designs, &actions, 1e-4),
        recommend(&designs, &actions, 1e-4)
    );
}

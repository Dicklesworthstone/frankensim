//! Battery for objective epistemics (addendum Proposal F). Covers CVaR, the
//! weakest-input color rule, robust-vs-nominal divergence, the amended
//! optimization contract (no optimizing an un-colored objective), the
//! kill-criterion dominance test, and colored fragility curves.

use fs_robust::{
    Color, ColorRank, ColoredObjective, RobustError, cvar, dominated_by_nominal, fragility_curve,
    robust_optimum, weakest_color,
};

fn verified() -> Color {
    Color::Verified { lo: -1.0, hi: 1.0 }
}
fn estimated() -> Color {
    Color::Estimated {
        estimator: "hazard-surrogate".into(),
        dispersion: 5.0,
    }
}

#[test]
fn cvar_weights_the_worst_tail() {
    let samples: Vec<f64> = (1..=100).map(f64::from).collect();
    // worst 10% (91..=100) has mean 95.5.
    assert!((cvar(&samples, 0.9).unwrap() - 95.5).abs() < 1e-9);
    // worst 5% (96..=100) has mean 98.0.
    assert!((cvar(&samples, 0.95).unwrap() - 98.0).abs() < 1e-9);
    // CVaR of the tail is worse than the mean (50.5).
    assert!(cvar(&samples, 0.9).unwrap() > 50.5);
}

#[test]
fn cvar_fractionally_weights_a_non_integral_tail_boundary() {
    // n*(1-alpha) = 1.5: the worst sample has full mass and the boundary
    // sample has half mass. Equal-weighting the top two would return 15 and
    // under-report this upper-tail risk.
    let actual = cvar(&[0.0, 10.0, 20.0], 0.5).unwrap();
    let expected = (20.0 + 0.5 * 10.0) / 1.5;
    assert!((actual - expected).abs() < 1e-12, "{actual} vs {expected}");
    assert!(
        actual > 15.0,
        "the old rounded-tail estimator was anti-conservative"
    );
}

#[test]
fn risk_means_do_not_overflow_on_finite_constant_samples() {
    let samples = [f64::MAX, f64::MAX, f64::MAX];
    assert_eq!(cvar(&samples, 0.5).unwrap().to_bits(), f64::MAX.to_bits());
    let objective = ColoredObjective::new("extreme", samples.to_vec(), vec![verified()]);
    assert_eq!(
        objective.nominal_value().unwrap().to_bits(),
        f64::MAX.to_bits()
    );
    let mixed = ColoredObjective::new(
        "mixed-extremes",
        vec![f64::MAX, -f64::MAX],
        vec![verified()],
    );
    let reversed = ColoredObjective::new(
        "mixed-extremes-reversed",
        vec![-f64::MAX, f64::MAX],
        vec![verified()],
    );
    assert_eq!(mixed.nominal_value().unwrap().to_bits(), 0.0_f64.to_bits());
    assert_eq!(
        mixed.nominal_value().unwrap().to_bits(),
        reversed.nominal_value().unwrap().to_bits(),
        "a sample statistic must not depend on input permutation"
    );
}

#[test]
fn cvar_rejects_bad_inputs() {
    assert_eq!(cvar(&[], 0.9), Err(RobustError::EmptySamples));
    assert!(matches!(
        cvar(&[1.0], 0.0),
        Err(RobustError::BadAlpha { .. })
    ));
    assert!(matches!(
        cvar(&[1.0], 1.0),
        Err(RobustError::BadAlpha { .. })
    ));
    assert!(matches!(
        cvar(&[1.0], 1.5),
        Err(RobustError::BadAlpha { .. })
    ));
    assert!(matches!(
        cvar(&[1.0], f64::NAN),
        Err(RobustError::BadAlpha { .. })
    ));
    assert!(matches!(
        cvar(&[1.0, f64::INFINITY], 0.9),
        Err(RobustError::BadSample { value }) if value.is_infinite()
    ));
    assert!(matches!(
        cvar(&[1.0, f64::NAN], 0.9),
        Err(RobustError::BadSample { value }) if value.is_nan()
    ));
}

#[test]
fn the_headline_takes_the_weakest_input_color() {
    // a verified structural solve under an estimated hazard is ESTIMATED.
    assert_eq!(
        weakest_color(&[verified(), estimated()]).unwrap().rank(),
        ColorRank::Estimated
    );
    let obj = ColoredObjective::new("d", vec![1.0, 2.0, 3.0], vec![verified(), estimated()]);
    assert_eq!(obj.headline_color().unwrap().rank(), ColorRank::Estimated);
    assert!(weakest_color(&[]).is_none());
}

#[test]
fn robust_and_nominal_optima_can_diverge() {
    // design A: low mean (10.8) but a fat tail (CVaR 50).
    let a = ColoredObjective::new("A", vec![1.0, 1.0, 1.0, 1.0, 50.0], vec![verified()]);
    // design B: higher mean (12) but tight (CVaR 12).
    let b = ColoredObjective::new("B", vec![12.0, 12.0, 12.0, 12.0, 12.0], vec![verified()]);
    // nominal (mean) would prefer A...
    assert!(a.nominal_value().unwrap() < b.nominal_value().unwrap());
    // ...but the ROBUST optimum (min CVaR) prefers B.
    let report = robust_optimum(&[a, b], 0.8).unwrap();
    assert_eq!(report.design, "B");
    assert!((report.robust_value - 12.0).abs() < 1e-9);
}

#[test]
fn optimization_refuses_an_un_colored_objective() {
    // the amended optimization contract: no color -> no optimization.
    let uncolored = ColoredObjective::new("fiction", vec![1.0, 2.0], vec![]);
    assert!(matches!(
        robust_optimum(std::slice::from_ref(&uncolored), 0.9),
        Err(RobustError::UncoloredObjective { design }) if design == "fiction"
    ));
    assert!(matches!(
        uncolored.headline_color(),
        Err(RobustError::UncoloredObjective { .. })
    ));
    // no candidates at all is also refused.
    assert_eq!(robust_optimum(&[], 0.9), Err(RobustError::NoCandidates));
    let bad_samples = ColoredObjective::new("bad", vec![1.0, f64::NAN], vec![verified()]);
    assert!(matches!(
        bad_samples.nominal_value(),
        Err(RobustError::BadSample { value }) if value.is_nan()
    ));
    assert!(matches!(
        robust_optimum(&[bad_samples], 0.9),
        Err(RobustError::BadSample { value }) if value.is_nan()
    ));
}

#[test]
fn the_kill_criterion_detects_domination_by_nominal_plus_safety() {
    // robust design costs 100, nominal+safety costs 90 -> robust is dominated.
    assert!(dominated_by_nominal(100.0, 90.0).unwrap());
    // robust costs 80, nominal+safety costs 90 -> robust wins (not dominated).
    assert!(!dominated_by_nominal(80.0, 90.0).unwrap());
    assert!(matches!(
        dominated_by_nominal(f64::NAN, 90.0),
        Err(RobustError::BadSample { value }) if value.is_nan()
    ));
    assert!(matches!(
        dominated_by_nominal(100.0, f64::INFINITY),
        Err(RobustError::BadSample { value }) if value.is_infinite()
    ));
}

#[test]
fn fragility_curves_are_monotone_and_colored() {
    // capacities clustered near 5; failure = demand exceeds capacity.
    let capacities = vec![3.0, 4.0, 5.0, 6.0, 7.0];
    let intensities = vec![1.0, 4.0, 6.0, 9.0];
    let f = fragility_curve(&capacities, &intensities, estimated()).unwrap();
    // P(failure) is nondecreasing in intensity, 0 at low, 1 at high.
    assert!((f.curve[0].prob_failure - 0.0).abs() < 1e-12);
    assert!((f.curve[3].prob_failure - 1.0).abs() < 1e-12);
    for w in f.curve.windows(2) {
        assert!(w[1].prob_failure >= w[0].prob_failure);
    }
    // the curve carries its honest color band.
    assert_eq!(f.color.rank(), ColorRank::Estimated);
    assert_eq!(
        fragility_curve(&[], &intensities, verified()),
        Err(RobustError::EmptySamples)
    );
    assert!(matches!(
        fragility_curve(&[3.0, f64::NAN], &intensities, verified()),
        Err(RobustError::BadSample { value }) if value.is_nan()
    ));
    assert!(matches!(
        fragility_curve(&capacities, &[1.0, f64::INFINITY], verified()),
        Err(RobustError::BadSample { value }) if value.is_infinite()
    ));
}

#[test]
fn fragility_curve_canonicalizes_unsorted_intensities() {
    let capacities = [3.0, 4.0, 5.0, 6.0, 7.0];
    let result = fragility_curve(&capacities, &[9.0, 1.0, 6.0, 4.0], estimated()).unwrap();
    let intensities: Vec<f64> = result.curve.iter().map(|point| point.intensity).collect();
    assert_eq!(intensities, vec![1.0, 4.0, 6.0, 9.0]);
    assert!(
        result
            .curve
            .windows(2)
            .all(|pair| pair[0].prob_failure <= pair[1].prob_failure)
    );
}

#[test]
fn optimization_is_deterministic() {
    let a = ColoredObjective::new("A", vec![1.0, 1.0, 50.0], vec![verified()]);
    let b = ColoredObjective::new("B", vec![12.0, 12.0, 12.0], vec![verified()]);
    assert_eq!(
        robust_optimum(&[a.clone(), b.clone()], 0.8),
        robust_optimum(&[a, b], 0.8)
    );
}

// ── Admitted-headline battery (bead 6pf9, stage S2) ─────────────────────────

use fs_evidence::{
    AdmissionDecision, AdmissionReceipt, AdmissionVerifier, AdmittedColor, COLOR_ALGEBRA_VERSION,
};
use fs_robust::{AdmittedRobustReport, admitted_headline_for, robust_optimum_admitted};

/// Test-fixture capability: accepts everything, standing in for the ledger
/// authority at the composition root. The real deny-all default and the
/// ledger oracle are exercised in fs-evidence and fs-ledger batteries.
struct FixtureAuthority;
impl AdmissionVerifier for FixtureAuthority {
    fn verify(&self, _c: &Color, _r: &AdmissionReceipt) -> AdmissionDecision {
        AdmissionDecision::accept(fs_blake3::hash_bytes(b"fixture-authority"))
    }
}

fn admit(color: Color, tag: &[u8]) -> AdmittedColor {
    AdmittedColor::from_receipt(
        color,
        AdmissionReceipt::from_parts(
            fs_blake3::hash_bytes(tag),
            7,
            COLOR_ALGEBRA_VERSION,
            fs_blake3::hash_bytes(b"fixture-policy"),
        ),
        &FixtureAuthority,
    )
    .expect("fixture admission")
}

fn verified_wide() -> Color {
    Color::Verified { lo: -9.0, hi: 9.0 }
}

#[test]
fn weakest_color_is_permutation_invariant_on_rank_ties() {
    let forward = weakest_color(&[verified(), verified_wide()]).expect("tie");
    let reversed = weakest_color(&[verified_wide(), verified()]).expect("tie");
    assert_eq!(
        forward, reversed,
        "reordering equal-rank inputs must not change the reported payload"
    );
    // The winner is content-determined (canonical-bytes minimum), so a
    // shuffled triple agrees too.
    let triple = weakest_color(&[verified_wide(), verified(), verified_wide()]).expect("tie");
    assert_eq!(triple, forward);
}

#[test]
fn headline_refuses_structurally_malformed_inputs() {
    let obj = ColoredObjective::new(
        "garbage",
        vec![1.0, 2.0],
        vec![Color::Verified {
            lo: f64::NAN,
            hi: 1.0,
        }],
    );
    assert!(matches!(
        obj.headline_color(),
        Err(RobustError::MalformedInputColor { design, .. }) if design == "garbage"
    ));
}

#[test]
fn admitted_headline_requires_full_count_aware_coverage() {
    let objective = ColoredObjective::new(
        "bridge",
        vec![1.0, 2.0],
        vec![verified(), verified(), verified_wide()],
    );
    // Full coverage: two admitted copies of the duplicate + the wide one.
    let full = [
        admit(verified(), b"node-a"),
        admit(verified(), b"node-b"),
        admit(verified_wide(), b"node-c"),
    ];
    let headline = admitted_headline_for(&objective, &full).expect("covered");
    assert_eq!(headline.rank(), ColorRank::Verified);

    // Count-aware: ONE admitted copy cannot cover TWO identical declared
    // inputs.
    let short = [
        admit(verified(), b"node-a"),
        admit(verified_wide(), b"node-c"),
    ];
    assert!(matches!(
        admitted_headline_for(&objective, &short),
        Err(RobustError::UnadmittedInput { design }) if design == "bridge"
    ));

    // An estimated declared input can never be covered: AdmittedColor is
    // always positive, so the objective keeps a declared-only headline.
    let mixed = ColoredObjective::new("mixed", vec![1.0], vec![verified(), estimated()]);
    assert!(matches!(
        admitted_headline_for(&mixed, &full),
        Err(RobustError::UnadmittedInput { design }) if design == "mixed"
    ));
}

#[test]
fn admitted_headline_ignores_surplus_admitted_values_and_is_order_free() {
    // The caller holds an admitted color WEAKER-ranked than anything the
    // objective declares; it must not leak into the headline.
    let regime = fs_evidence::ValidityDomain::unconstrained().with("re", 1e3, 1e5);
    let surplus = admit(
        Color::Validated {
            regime,
            dataset: "unrelated-campaign".to_string(),
        },
        b"node-surplus",
    );
    let objective = ColoredObjective::new("clean", vec![1.0], vec![verified()]);
    let forward = [admit(verified(), b"node-a"), surplus.clone()];
    let reversed = [surplus, admit(verified(), b"node-a")];
    let a = admitted_headline_for(&objective, &forward).expect("covered");
    let b = admitted_headline_for(&objective, &reversed).expect("covered");
    assert_eq!(a.rank(), ColorRank::Verified, "surplus Validated must not leak");
    assert_eq!(a, b, "admitted headline must be order-free");
}

#[test]
fn admitted_headline_same_color_surplus_is_order_free() {
    // Regression: a SAME-COLOR surplus — two admitted counterparts that both
    // cover the single declared input but carry DIFFERENT receipts. Only one is
    // consumed (count-aware coverage), and WHICH one must be canonical, not
    // input-order-dependent. The old first-match consumption let the headline
    // lineage flip on reordering; the existing surplus test missed it because
    // its surplus was a DIFFERENT color (never a coverage candidate).
    let objective = ColoredObjective::new("clean", vec![1.0], vec![verified()]);
    let node_a = admit(verified(), b"node-a");
    let node_b = admit(verified(), b"node-b");
    let forward = [node_a.clone(), node_b.clone()];
    let reversed = [node_b, node_a];
    let a = admitted_headline_for(&objective, &forward).expect("covered");
    let b = admitted_headline_for(&objective, &reversed).expect("covered");
    assert_eq!(
        a, b,
        "same-color surplus: the headline receipt must not depend on input order"
    );
}

#[test]
fn admitted_optimum_requires_every_candidate_admitted() {
    let a = ColoredObjective::new("A", vec![1.0, 1.0, 1.0, 1.0, 50.0], vec![verified()]);
    let b = ColoredObjective::new("B", vec![12.0; 5], vec![verified()]);
    let admitted_a = [admit(verified(), b"node-a")];
    let admitted_b = [admit(verified(), b"node-b")];

    let report: AdmittedRobustReport = robust_optimum_admitted(
        &[
            (a.clone(), admitted_a.as_slice()),
            (b.clone(), admitted_b.as_slice()),
        ],
        0.8,
    )
    .expect("wholly admitted run");
    assert_eq!(report.design, "B");
    assert_eq!(report.headline.rank(), ColorRank::Verified);
    assert_eq!(
        report.headline.receipt().node_hash(),
        fs_blake3::hash_bytes(b"node-b"),
        "the headline must carry the WINNER'S admission lineage"
    );

    // One unadmitted candidate poisons the whole positive claim — even when
    // that candidate would lose the optimization anyway.
    assert!(matches!(
        robust_optimum_admitted(&[(a, admitted_a.as_slice()), (b, &[])], 0.8),
        Err(RobustError::UnadmittedInput { design }) if design == "B"
    ));
}

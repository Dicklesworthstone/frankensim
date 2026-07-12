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

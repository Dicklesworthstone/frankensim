//! Battery for adjoint-driven tolerance allocation (addendum Proposal 11).
//! Covers the tighten-high / loosen-low allocation meeting the variance budget,
//! the band-extremes robustness check, the GD&T report carrying certified
//! sensitivities, the P(in-spec) → variance budget, and error paths.

use fs_toleralloc::{
    Action, ColorRank, Feature, ToleranceError, allocate, gdt_report, robustness_check,
    variance_budget,
};

fn feature(name: &str, sensitivity: f64, baseline: f64) -> Feature {
    Feature {
        name: name.into(),
        sensitivity,
        sensitivity_color: ColorRank::Verified,
        cost_coeff: 1.0,
        baseline_tolerance: baseline,
    }
}

#[test]
fn tolerance_is_spent_where_sensitivity_is_large() {
    // a critical (high-sensitivity) feature and a slack (low-sensitivity) one.
    let features = vec![feature("critical", 10.0, 0.5), feature("slack", 0.1, 0.5)];
    let alloc = allocate(&features, 1.0, 3.0).unwrap();
    let crit = &alloc.items[0];
    let slack = &alloc.items[1];
    // the high-sensitivity feature is tightened; the low-sensitivity one loosened.
    assert_eq!(crit.action, Action::Tighten, "crit tol {}", crit.tolerance);
    assert_eq!(
        slack.action,
        Action::Loosen,
        "slack tol {}",
        slack.tolerance
    );
    assert!(crit.tolerance < slack.tolerance);
    // the budget is met exactly (by construction).
    assert!((alloc.achieved_variance - 1.0).abs() < 1e-9);
}

#[test]
fn allocation_rejects_bad_input() {
    assert_eq!(allocate(&[], 1.0, 3.0), Err(ToleranceError::NoFeatures));
    assert!(matches!(
        allocate(&[feature("z", 0.0, 0.5)], 1.0, 3.0),
        Err(ToleranceError::NonPositive {
            what: "sensitivity",
            ..
        })
    ));
    assert_eq!(
        allocate(&[feature("f", 1.0, 0.5)], 0.0, 3.0),
        Err(ToleranceError::BadBudget)
    );
    assert_eq!(
        allocate(&[feature("f", 1.0, 0.5)], 1.0, 0.0),
        Err(ToleranceError::BadBudget)
    );
}

#[test]
fn the_robustness_check_confirms_or_flags_the_linearization() {
    let alloc = allocate(&[feature("a", 1.0, 0.5), feature("b", 2.0, 0.5)], 1.0, 3.0).unwrap();
    // linearized std = sqrt(budget) = 1; bound = 3 * 1 * 1.2 = 3.6.
    // extremes within the bound -> confirmed.
    let ok = robustness_check(&alloc, &[2.5, -2.0, 1.0], 0.0, 3.0, 0.2);
    assert!(ok.confirmed);
    assert!((ok.linearized_std - 1.0).abs() < 1e-9);
    // an extreme far beyond the linear prediction -> flagged (nonlinearity).
    let bad = robustness_check(&alloc, &[8.0], 0.0, 3.0, 0.2);
    assert!(!bad.confirmed);
    assert!((bad.sampled_max_deviation - 8.0).abs() < 1e-12);
}

#[test]
fn the_gdt_report_attaches_a_certified_sensitivity_to_every_loosened_tolerance() {
    let features = vec![feature("critical", 10.0, 0.5), feature("slack", 0.1, 0.5)];
    let alloc = allocate(&features, 1.0, 3.0).unwrap();
    let report = gdt_report(&alloc);
    assert_eq!(report.len(), 2);
    for s in &report {
        // every suggestion carries the certified sensitivity + its color.
        assert!(s.certified_sensitivity > 0.0);
        assert_eq!(s.color, ColorRank::Verified);
    }
    // the loosened tolerance (the savings) is justified by its low sensitivity.
    let loosened: Vec<_> = report
        .iter()
        .filter(|s| s.action == Action::Loosen)
        .collect();
    assert_eq!(loosened.len(), 1);
    assert!((loosened[0].certified_sensitivity - 0.1).abs() < 1e-12);
}

#[test]
fn the_variance_budget_follows_the_in_spec_probability() {
    // P(|QoI| <= 1.96 sigma) = 0.95 => sigma = 1 => budget = 1.
    let b = variance_budget(1.96, 0.95).unwrap();
    assert!((b - 1.0).abs() < 1e-2, "budget {b}");
    // a tighter target needs a smaller budget (smaller allowed variance).
    let tight = variance_budget(1.0, 0.99).unwrap();
    let loose = variance_budget(1.0, 0.90).unwrap();
    assert!(tight < loose);
    // bad inputs.
    assert_eq!(variance_budget(1.0, 1.0), Err(ToleranceError::BadBudget));
    assert_eq!(variance_budget(0.0, 0.95), Err(ToleranceError::BadBudget));
}

#[test]
fn allocation_is_deterministic() {
    let features = vec![
        feature("a", 3.0, 0.2),
        feature("b", 0.5, 0.2),
        feature("c", 1.0, 0.2),
    ];
    assert_eq!(allocate(&features, 0.5, 3.0), allocate(&features, 0.5, 3.0));
}

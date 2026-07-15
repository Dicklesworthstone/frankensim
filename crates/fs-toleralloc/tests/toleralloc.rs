//! Battery for adjoint-driven tolerance allocation (addendum Proposal 11).
//! Covers the tighten-high / loosen-low allocation meeting the variance budget,
//! the band-extremes robustness check, the GD&T report carrying certified
//! sensitivities, the P(in-spec) → variance budget, and error paths.

use fs_toleralloc::{
    Action, Allocation, ColorRank, DerivedQuantity, Feature, ScalarIssue, ToleranceError, allocate,
    gdt_report, robustness_check, variance_budget,
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
        Err(ToleranceError::InvalidFeatureField {
            index: 0,
            field: "sensitivity",
            issue: ScalarIssue::NonPositive,
            ..
        })
    ));
    assert!(matches!(
        allocate(&[feature("f", 1.0, 0.5)], 0.0, 3.0),
        Err(ToleranceError::InvalidArgument {
            argument: "variance_budget",
            issue: ScalarIssue::NonPositive,
        })
    ));
    assert!(matches!(
        allocate(&[feature("f", 1.0, 0.5)], 1.0, 0.0),
        Err(ToleranceError::InvalidArgument {
            argument: "k",
            issue: ScalarIssue::NonPositive,
        })
    ));

    for value in [0.0, -1.0] {
        let mut bad = feature("s", value, 0.5);
        assert!(matches!(
            allocate(&[bad.clone()], 1.0, 3.0),
            Err(ToleranceError::InvalidFeatureField {
                field: "sensitivity",
                issue: ScalarIssue::NonPositive,
                ..
            })
        ));
        bad.sensitivity = 1.0;
        bad.cost_coeff = value;
        assert!(matches!(
            allocate(&[bad.clone()], 1.0, 3.0),
            Err(ToleranceError::InvalidFeatureField {
                field: "cost_coeff",
                issue: ScalarIssue::NonPositive,
                ..
            })
        ));
        bad.cost_coeff = 1.0;
        bad.baseline_tolerance = value;
        assert!(matches!(
            allocate(&[bad], 1.0, 3.0),
            Err(ToleranceError::InvalidFeatureField {
                field: "baseline_tolerance",
                issue: ScalarIssue::NonPositive,
                ..
            })
        ));
        assert!(matches!(
            allocate(&[feature("f", 1.0, 0.5)], value, 3.0),
            Err(ToleranceError::InvalidArgument {
                argument: "variance_budget",
                issue: ScalarIssue::NonPositive,
            })
        ));
        assert!(matches!(
            allocate(&[feature("f", 1.0, 0.5)], 1.0, value),
            Err(ToleranceError::InvalidArgument {
                argument: "k",
                issue: ScalarIssue::NonPositive,
            })
        ));
    }
}

#[test]
fn the_robustness_check_confirms_or_flags_the_linearization() {
    let alloc = allocate(&[feature("a", 1.0, 0.5), feature("b", 2.0, 0.5)], 1.0, 3.0).unwrap();
    // linearized std = sqrt(budget) = 1; bound = 3 * 1 * 1.2 = 3.6.
    // extremes within the bound -> confirmed.
    let ok = robustness_check(&alloc, &[2.5, -2.0, 1.0], 0.0, 3.0, 0.2).unwrap();
    assert!(ok.confirmed);
    assert!((ok.linearized_std - 1.0).abs() < 1e-9);
    // an extreme far beyond the linear prediction -> flagged (nonlinearity).
    let bad = robustness_check(&alloc, &[8.0], 0.0, 3.0, 0.2).unwrap();
    assert!(!bad.confirmed);
    assert!((bad.sampled_max_deviation - 8.0).abs() < 1e-12);
}

#[test]
fn the_gdt_report_attaches_a_certified_sensitivity_to_every_loosened_tolerance() {
    let features = vec![feature("critical", 10.0, 0.5), feature("slack", 0.1, 0.5)];
    let alloc = allocate(&features, 1.0, 3.0).unwrap();
    let report = gdt_report(&alloc).unwrap();
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
    assert!(matches!(
        variance_budget(1.0, 1.0),
        Err(ToleranceError::InvalidArgument {
            argument: "target",
            issue: ScalarIssue::OutsideOpenUnitInterval,
        })
    ));
    assert!(matches!(
        variance_budget(0.0, 0.95),
        Err(ToleranceError::InvalidArgument {
            argument: "spec_margin",
            issue: ScalarIssue::NonPositive,
        })
    ));

    // Adjacent representable targets must not round the internal CDF input to
    // exactly 0.5 or 1.0 before evaluating the quantile.
    let near_one = variance_budget(1.0, 1.0_f64.next_down()).unwrap();
    assert!(near_one.is_finite() && near_one > 0.0);
    let near_zero = variance_budget(f64::MIN_POSITIVE, f64::from_bits(1)).unwrap();
    assert!(near_zero.is_finite() && near_zero > 0.0);
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

#[test]
fn every_non_finite_public_input_is_refused_at_its_field() {
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut bad = feature("s", 1.0, 0.5);
        bad.sensitivity = value;
        assert!(matches!(
            allocate(&[bad], 1.0, 3.0),
            Err(ToleranceError::InvalidFeatureField {
                index: 0,
                field: "sensitivity",
                issue: ScalarIssue::NonFinite,
                ..
            })
        ));

        let mut bad = feature("c", 1.0, 0.5);
        bad.cost_coeff = value;
        assert!(matches!(
            allocate(&[bad], 1.0, 3.0),
            Err(ToleranceError::InvalidFeatureField {
                index: 0,
                field: "cost_coeff",
                issue: ScalarIssue::NonFinite,
                ..
            })
        ));

        let mut bad = feature("b", 1.0, 0.5);
        bad.baseline_tolerance = value;
        assert!(matches!(
            allocate(&[bad], 1.0, 3.0),
            Err(ToleranceError::InvalidFeatureField {
                index: 0,
                field: "baseline_tolerance",
                issue: ScalarIssue::NonFinite,
                ..
            })
        ));

        assert!(matches!(
            allocate(&[feature("f", 1.0, 0.5)], value, 3.0),
            Err(ToleranceError::InvalidArgument {
                argument: "variance_budget",
                issue: ScalarIssue::NonFinite,
            })
        ));
        assert!(matches!(
            allocate(&[feature("f", 1.0, 0.5)], 1.0, value),
            Err(ToleranceError::InvalidArgument {
                argument: "k",
                issue: ScalarIssue::NonFinite,
            })
        ));
        assert!(matches!(
            variance_budget(value, 0.95),
            Err(ToleranceError::InvalidArgument {
                argument: "spec_margin",
                issue: ScalarIssue::NonFinite,
            })
        ));
        assert!(matches!(
            variance_budget(1.0, value),
            Err(ToleranceError::InvalidArgument {
                argument: "target",
                issue: ScalarIssue::NonFinite,
            })
        ));
    }
}

#[test]
fn feature_names_are_stable_and_unambiguous() {
    assert!(matches!(
        allocate(&[feature("", 1.0, 0.5)], 1.0, 3.0),
        Err(ToleranceError::InvalidFeatureName {
            index: 0,
            reason: "name must not be empty",
            ..
        })
    ));
    assert!(matches!(
        allocate(&[feature(" edge", 1.0, 0.5)], 1.0, 3.0),
        Err(ToleranceError::InvalidFeatureName {
            index: 0,
            reason: "name must not have leading or trailing whitespace",
            ..
        })
    ));
    assert!(matches!(
        allocate(&[feature("edge\u{0000}face", 1.0, 0.5)], 1.0, 3.0),
        Err(ToleranceError::InvalidFeatureName {
            index: 0,
            reason: "name must not contain control characters",
            ..
        })
    ));
    assert!(matches!(
        allocate(
            &[feature("edge", 1.0, 0.5), feature("edge", 2.0, 0.5)],
            1.0,
            3.0,
        ),
        Err(ToleranceError::AmbiguousFeatureName {
            first_index: 0,
            duplicate_index: 1,
            ref canonical_name,
        }) if canonical_name == "edge"
    ));
    assert!(matches!(
        allocate(
            &[feature("LeadingEdge", 1.0, 0.5), feature("leadingedge", 2.0, 0.5)],
            1.0,
            3.0,
        ),
        Err(ToleranceError::AmbiguousFeatureName {
            first_index: 0,
            duplicate_index: 1,
            ref canonical_name,
        }) if canonical_name == "leadingedge"
    ));
}

#[test]
fn boundary_values_never_publish_non_finite_outputs() {
    let tiny = Feature {
        name: "tiny".into(),
        sensitivity: f64::MIN_POSITIVE,
        sensitivity_color: ColorRank::Verified,
        cost_coeff: f64::MIN_POSITIVE,
        baseline_tolerance: f64::MIN_POSITIVE,
    };
    let allocation = allocate(&[tiny], f64::MIN_POSITIVE, f64::MIN_POSITIVE).unwrap();
    assert!(allocation.total_cost.is_finite() && allocation.total_cost > 0.0);
    assert!(allocation.achieved_variance.is_finite() && allocation.achieved_variance > 0.0);
    assert!(
        allocation
            .items
            .iter()
            .all(|item| item.tolerance.is_finite() && item.tolerance > 0.0)
    );

    let huge_cost = Feature {
        name: "unrepresentable-cost".into(),
        sensitivity: f64::MAX,
        sensitivity_color: ColorRank::Verified,
        cost_coeff: f64::MAX,
        baseline_tolerance: 1.0,
    };
    assert!(matches!(
        allocate(&[huge_cost], 1.0, 1.0),
        Err(ToleranceError::InvalidDerived {
            quantity: DerivedQuantity::CostContribution,
            feature_index: Some(0),
            issue: ScalarIssue::NonFinite,
        })
    ));
}

#[test]
fn robustness_refuses_empty_poisoned_and_unrepresentable_evidence() {
    let alloc = allocate(&[feature("a", 1.0, 0.5)], 1.0, 3.0).unwrap();
    assert_eq!(
        robustness_check(&alloc, &[], 0.0, 3.0, 0.0),
        Err(ToleranceError::NoExtremeSamples)
    );
    assert!(matches!(
        robustness_check(&alloc, &[f64::NAN], 0.0, 3.0, 0.0),
        Err(ToleranceError::InvalidExtremeQoi {
            index: 0,
            issue: ScalarIssue::NonFinite,
        })
    ));
    assert!(matches!(
        robustness_check(&alloc, &[0.0], f64::INFINITY, 3.0, 0.0),
        Err(ToleranceError::InvalidArgument {
            argument: "nominal_qoi",
            issue: ScalarIssue::NonFinite,
        })
    ));
    assert!(matches!(
        robustness_check(&alloc, &[0.0], 0.0, 3.0, -0.1),
        Err(ToleranceError::InvalidArgument {
            argument: "margin",
            issue: ScalarIssue::Negative,
        })
    ));
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            robustness_check(&alloc, &[0.0], 0.0, value, 0.0),
            Err(ToleranceError::InvalidArgument {
                argument: "k",
                issue: ScalarIssue::NonFinite,
            })
        ));
        assert!(matches!(
            robustness_check(&alloc, &[0.0], 0.0, 3.0, value),
            Err(ToleranceError::InvalidArgument {
                argument: "margin",
                issue: ScalarIssue::NonFinite,
            })
        ));
    }

    let mut poisoned = alloc.clone();
    poisoned.achieved_variance = f64::NAN;
    assert!(matches!(
        robustness_check(&poisoned, &[0.0], 0.0, 3.0, 0.0),
        Err(ToleranceError::InvalidArgument {
            argument: "allocation.achieved_variance",
            issue: ScalarIssue::NonFinite,
        })
    ));

    assert!(matches!(
        robustness_check(&alloc, &[f64::MAX], -f64::MAX, 3.0, 0.0),
        Err(ToleranceError::InvalidExtremeDerived {
            index: 0,
            quantity: DerivedQuantity::SampledDeviation,
            issue: ScalarIssue::NonFinite,
        })
    ));
}

#[test]
fn gdt_report_refuses_forged_non_finite_items() {
    let mut allocation = allocate(&[feature("edge", 1.0, 0.5)], 1.0, 3.0).unwrap();
    allocation.items[0].tolerance = f64::NAN;
    assert!(matches!(
        gdt_report(&allocation),
        Err(ToleranceError::InvalidAllocationItem {
            index: 0,
            field: "tolerance",
            issue: ScalarIssue::NonFinite,
            ..
        })
    ));

    let mut allocation = allocate(&[feature("edge", 1.0, 0.5)], 1.0, 3.0).unwrap();
    allocation.items[0].sensitivity = f64::INFINITY;
    assert!(matches!(
        gdt_report(&allocation),
        Err(ToleranceError::InvalidAllocationItem {
            index: 0,
            field: "sensitivity",
            issue: ScalarIssue::NonFinite,
            ..
        })
    ));
}

#[test]
fn g3_common_sensitivity_rescaling_preserves_tolerances() {
    let original = vec![feature("a", 2.0, 0.5), feature("b", 5.0, 0.5)];
    let rescaled = vec![feature("a", 20.0, 0.5), feature("b", 50.0, 0.5)];
    let base = allocate(&original, 0.25, 3.0).unwrap();
    let scaled = allocate(&rescaled, 25.0, 3.0).unwrap();
    for (left, right) in base.items.iter().zip(&scaled.items) {
        let relative = (left.tolerance - right.tolerance).abs() / left.tolerance;
        assert!(
            relative < 1e-12,
            "{} rescaling drift: {relative}",
            left.name
        );
        assert_eq!(left.action, right.action);
    }
}

#[test]
fn g5_input_order_is_the_stable_tie_break() {
    let tied = vec![
        feature("zeta", 1.0, 0.5),
        feature("alpha", 1.0, 0.5),
        feature("middle", 1.0, 0.5),
    ];
    let first = allocate(&tied, 0.5, 3.0).unwrap();
    let second = allocate(&tied, 0.5, 3.0).unwrap();
    assert_eq!(first, second);
    assert_eq!(
        first
            .items
            .iter()
            .map(|item| item.name.as_str())
            .collect::<Vec<_>>(),
        ["zeta", "alpha", "middle"]
    );
}

#[test]
fn forged_empty_or_zero_allocation_cannot_confirm_or_publish() {
    let allocation = Allocation {
        items: Vec::new(),
        total_cost: 0.0,
        achieved_variance: 0.0,
    };
    assert_eq!(
        robustness_check(&allocation, &[0.0], 0.0, 1.0, 0.0),
        Err(ToleranceError::NoFeatures)
    );
    assert_eq!(gdt_report(&allocation), Err(ToleranceError::NoFeatures));

    let mut allocation = allocate(&[feature("edge", 1.0, 0.5)], 1.0, 3.0).unwrap();
    allocation.achieved_variance = 0.0;
    assert!(matches!(
        robustness_check(&allocation, &[0.0], 0.0, 1.0, 0.0),
        Err(ToleranceError::InvalidArgument {
            argument: "allocation.achieved_variance",
            issue: ScalarIssue::NonPositive,
        })
    ));
}

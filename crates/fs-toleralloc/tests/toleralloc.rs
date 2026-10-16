//! Battery for adjoint-driven tolerance allocation (addendum Proposal 11).
//! Covers the tighten-high / loosen-low allocation meeting the variance budget,
//! the band-extremes robustness check, the GD&T report carrying certified
//! sensitivities, the P(in-spec) → variance budget, and error paths.

use std::num::NonZeroU64;

use fs_toleralloc::{
    Action, AdmittedCorrelationModel, Allocation, ColorRank, CorrelatedDerivedQuantity,
    CorrelatedStackError, CorrelatedStackTerm, CorrelationAdmissionError, CorrelationFactorIssue,
    DerivedQuantity, Feature, MAX_CORRELATED_STACK_TERM_NAME_BYTES_V1,
    MAX_CORRELATED_STACK_TERMS_V1, ScalarIssue, ToleranceError, allocate, gdt_report,
    propagate_correlated_stack, robustness_check, variance_budget,
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

fn correlation_model(rho: f64, residual: f64) -> AdmittedCorrelationModel {
    AdmittedCorrelationModel::try_new(
        "gear/process-runout",
        NonZeroU64::new(1).expect("one is nonzero"),
        [0x5a; 32],
        2,
        vec![1.0, 0.0, rho, residual],
    )
    .expect("manufactured factor is admissible")
}

fn stack_term(name: &str, signed_sensitivity: f64, standard_deviation: f64) -> CorrelatedStackTerm {
    CorrelatedStackTerm {
        name: name.into(),
        signed_sensitivity,
        sensitivity_color: ColorRank::Verified,
        standard_deviation,
    }
}

fn assert_relative_close(actual: f64, expected: f64, tolerance: f64) {
    let scale = actual.abs().max(expected.abs()).max(1.0);
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}

fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut value = *state;
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

fn correlated_rademacher_monte_carlo_variance(samples: u32) -> f64 {
    let mut state = 0x6765_6172_2d73_7461_u64;
    let mut sum = 0.0_f64;
    let mut sum_squares = 0.0_f64;
    for _ in 0..samples {
        let x = if splitmix64(&mut state) & 1 == 0 {
            -1.0
        } else {
            1.0
        };
        let y = if splitmix64(&mut state) % 10 == 0 {
            -x
        } else {
            x
        };
        let qoi = x + y;
        sum += qoi;
        sum_squares += qoi * qoi;
    }
    let count = f64::from(samples);
    let mean = sum / count;
    sum_squares / count - mean * mean
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

#[test]
fn correlated_stack_catches_independence_error_against_exhaustive_population() {
    let model = correlation_model(0.8, 0.6);
    let terms = [
        stack_term("carrier-runout", 1.0, 1.0),
        stack_term("gear-eccentricity", 1.0, 1.0),
    ];
    let receipt = propagate_correlated_stack(&model, &terms).expect("correlated stack evaluates");

    // Manufactured finite population: Y equals X in 18 of 20 equiprobable
    // outcomes and is -X in two, so Corr(X,Y)=0.8 exactly. Exhaustive
    // enumeration is a stronger oracle than a sampled Monte Carlo estimate.
    let weighted_outputs = [(2.0_f64, 9_u32), (-2.0, 9), (0.0, 1), (0.0, 1)];
    let population_variance = weighted_outputs
        .iter()
        .map(|(value, weight)| value * value * f64::from(*weight))
        .sum::<f64>()
        / 20.0;
    let monte_carlo_variance = correlated_rademacher_monte_carlo_variance(200_000);

    assert_relative_close(receipt.independent_variance(), 2.0, 1e-15);
    assert_relative_close(receipt.correlated_variance(), 3.6, 1e-15);
    assert_relative_close(receipt.correlated_variance(), population_variance, 1e-15);
    assert!((receipt.correlated_variance() - monte_carlo_variance).abs() < 0.02);
    assert!((receipt.independent_variance() - monte_carlo_variance).abs() > 1.0);
    assert_relative_close(receipt.correlation_variance_delta(), 1.6, 1e-15);
    assert_relative_close(
        receipt.independent_standard_deviation(),
        2.0_f64.sqrt(),
        1e-15,
    );
    assert_relative_close(
        receipt.correlated_standard_deviation(),
        3.6_f64.sqrt(),
        1e-15,
    );
    assert_eq!(receipt.model().namespace(), "gear/process-runout");
    assert_eq!(receipt.model().schema_version().get(), 1);
    assert_eq!(receipt.model().semantic_digest(), [0x5a; 32]);
    assert_eq!(receipt.model().lower_factor(), [1.0, 0.0, 0.8, 0.6]);
    assert_eq!(receipt.terms(), terms);
    assert_eq!(
        receipt
            .terms()
            .iter()
            .map(|term| term.name.as_str())
            .collect::<Vec<_>>(),
        ["carrier-runout", "gear-eccentricity"]
    );
    assert_eq!(
        propagate_correlated_stack(&model, &terms),
        propagate_correlated_stack(&model, &terms)
    );
    let rebound_model = AdmittedCorrelationModel::try_new(
        "gear/process-runout",
        NonZeroU64::new(1).expect("one is nonzero"),
        [0x5b; 32],
        2,
        vec![1.0, 0.0, 0.8, 0.6],
    )
    .expect("rebound model is structurally admissible");
    let rebound =
        propagate_correlated_stack(&rebound_model, &terms).expect("rebound stack evaluates");
    assert_eq!(rebound.correlated_variance(), receipt.correlated_variance());
    assert_ne!(rebound, receipt, "model identity is receipt-semantic");
}

#[test]
fn g3_signed_sensitivity_and_correlation_can_reduce_or_increase_variance() {
    let positive = correlation_model(0.8, 0.6);
    let opposite_terms = [
        stack_term("carrier-runout", 1.0, 1.0),
        stack_term("gear-eccentricity", -1.0, 1.0),
    ];
    let reduced =
        propagate_correlated_stack(&positive, &opposite_terms).expect("opposite signs evaluate");
    assert_relative_close(reduced.independent_variance(), 2.0, 1e-15);
    assert_relative_close(reduced.correlated_variance(), 0.4, 1e-15);
    assert_relative_close(reduced.correlation_variance_delta(), -1.6, 1e-15);

    let negative = correlation_model(-0.8, 0.6);
    let increased =
        propagate_correlated_stack(&negative, &opposite_terms).expect("negative rho evaluates");
    assert_relative_close(increased.independent_variance(), 2.0, 1e-15);
    assert_relative_close(increased.correlated_variance(), 3.6, 1e-15);
    assert_relative_close(increased.correlation_variance_delta(), 1.6, 1e-15);

    let perfectly_correlated = correlation_model(1.0, 0.0);
    assert_eq!(
        propagate_correlated_stack(&perfectly_correlated, &opposite_terms),
        Err(CorrelatedStackError::InvalidDerived {
            quantity: CorrelatedDerivedQuantity::CorrelationProjection,
            term_index: None,
            issue: ScalarIssue::AmbiguousZero,
        })
    );
}

#[test]
fn correlation_factor_admission_is_bounded_canonical_and_psd_by_construction() {
    let version = NonZeroU64::new(1).expect("one is nonzero");
    assert!(matches!(
        AdmittedCorrelationModel::try_new("Bad.Namespace", version, [1; 32], 1, vec![1.0]),
        Err(CorrelationAdmissionError::InvalidNamespace { .. })
    ));
    let overlong_namespace = "a".repeat(257);
    assert!(matches!(
        AdmittedCorrelationModel::try_new(
            overlong_namespace,
            version,
            [1; 32],
            1,
            vec![1.0],
        ),
        Err(CorrelationAdmissionError::InvalidNamespace {
            ref namespace,
            reason: "namespace exceeds the versioned byte cap",
        }) if namespace.len() == 256
    ));
    assert_eq!(
        AdmittedCorrelationModel::try_new("gear/runout", version, [0; 32], 1, vec![1.0]),
        Err(CorrelationAdmissionError::ZeroDigest)
    );
    assert_eq!(
        AdmittedCorrelationModel::try_new("gear/runout", version, [1; 32], 0, Vec::new()),
        Err(CorrelationAdmissionError::InvalidDimension {
            dimension: 0,
            max: MAX_CORRELATED_STACK_TERMS_V1,
        })
    );
    assert_eq!(
        AdmittedCorrelationModel::try_new(
            "gear/runout",
            version,
            [1; 32],
            MAX_CORRELATED_STACK_TERMS_V1 + 1,
            Vec::new(),
        ),
        Err(CorrelationAdmissionError::InvalidDimension {
            dimension: MAX_CORRELATED_STACK_TERMS_V1 + 1,
            max: MAX_CORRELATED_STACK_TERMS_V1,
        })
    );
    assert_eq!(
        AdmittedCorrelationModel::try_new("gear/runout", version, [1; 32], 2, vec![1.0; 3]),
        Err(CorrelationAdmissionError::FactorLength {
            dimension: 2,
            expected: 4,
            actual: 3,
        })
    );

    for (factor, row, column, issue) in [
        (
            vec![1.0, f64::NAN, 0.0, 1.0],
            0,
            1,
            CorrelationFactorIssue::NonFinite,
        ),
        (
            vec![1.0, -0.0, 0.0, 1.0],
            0,
            1,
            CorrelationFactorIssue::NonCanonicalNegativeZero,
        ),
        (
            vec![1.0, 0.25, 0.0, 1.0],
            0,
            1,
            CorrelationFactorIssue::AboveDiagonalNonZero,
        ),
        (
            vec![1.0, 0.0, 0.0, -1.0],
            1,
            1,
            CorrelationFactorIssue::NegativeDiagonal,
        ),
        (
            vec![1.0, 0.0, 1.01, 0.0],
            1,
            0,
            CorrelationFactorIssue::MagnitudeAboveOne,
        ),
    ] {
        assert_eq!(
            AdmittedCorrelationModel::try_new("gear/runout", version, [1; 32], 2, factor),
            Err(CorrelationAdmissionError::InvalidFactorEntry { row, column, issue })
        );
    }

    assert!(matches!(
        AdmittedCorrelationModel::try_new(
            "gear/runout",
            version,
            [1; 32],
            2,
            vec![1.0, 0.0, 0.5, 0.5],
        ),
        Err(CorrelationAdmissionError::NonUnitRow { row: 1, .. })
    ));
}

#[test]
fn correlated_stack_refuses_axis_and_numeric_ambiguity() {
    let model = correlation_model(0.8, 0.6);
    assert_eq!(
        propagate_correlated_stack(&model, &[]),
        Err(CorrelatedStackError::NoTerms)
    );
    assert_eq!(
        propagate_correlated_stack(&model, &[stack_term("only-one-axis", 1.0, 1.0)],),
        Err(CorrelatedStackError::DimensionMismatch { model: 2, terms: 1 })
    );
    let too_many = (0..=MAX_CORRELATED_STACK_TERMS_V1)
        .map(|index| stack_term(&format!("axis-{index}"), 1.0, 1.0))
        .collect::<Vec<_>>();
    assert_eq!(
        propagate_correlated_stack(&model, &too_many),
        Err(CorrelatedStackError::TooManyTerms {
            actual: MAX_CORRELATED_STACK_TERMS_V1 + 1,
            max: MAX_CORRELATED_STACK_TERMS_V1,
        })
    );
    assert!(matches!(
        propagate_correlated_stack(
            &model,
            &[
                stack_term("Runout", 1.0, 1.0),
                stack_term("runout", 1.0, 1.0),
            ],
        ),
        Err(CorrelatedStackError::AmbiguousTermName {
            first_index: 0,
            duplicate_index: 1,
            ..
        })
    ));
    assert!(matches!(
        propagate_correlated_stack(
            &model,
            &[
                stack_term(" bad-axis", 1.0, 1.0),
                stack_term("other", 1.0, 1.0),
            ],
        ),
        Err(CorrelatedStackError::InvalidTermName { index: 0, .. })
    ));
    let overlong_name = "a".repeat(MAX_CORRELATED_STACK_TERM_NAME_BYTES_V1 + 1);
    assert!(matches!(
        propagate_correlated_stack(
            &model,
            &[
                stack_term(&overlong_name, 1.0, 1.0),
                stack_term("other", 1.0, 1.0),
            ],
        ),
        Err(CorrelatedStackError::InvalidTermName {
            index: 0,
            ref name,
            reason: "name exceeds the versioned byte cap",
        }) if name.len() == MAX_CORRELATED_STACK_TERM_NAME_BYTES_V1
    ));
    let expanding_name = "\u{130}".repeat(MAX_CORRELATED_STACK_TERM_NAME_BYTES_V1 / 2);
    assert_eq!(
        expanding_name.len(),
        MAX_CORRELATED_STACK_TERM_NAME_BYTES_V1
    );
    assert!(matches!(
        propagate_correlated_stack(
            &model,
            &[
                stack_term(&expanding_name, 1.0, 1.0),
                stack_term("other", 1.0, 1.0),
            ],
        ),
        Err(CorrelatedStackError::InvalidTermName {
            index: 0,
            reason: "lowercase comparison key exceeds the versioned byte cap",
            ..
        })
    ));
    assert!(matches!(
        propagate_correlated_stack(
            &model,
            &[
                stack_term("bad-sensitivity", f64::INFINITY, 1.0),
                stack_term("other", 1.0, 1.0),
            ],
        ),
        Err(CorrelatedStackError::InvalidTermField {
            index: 0,
            field: "signed_sensitivity",
            issue: ScalarIssue::NonFinite,
            ..
        })
    ));
    assert!(matches!(
        propagate_correlated_stack(
            &model,
            &[
                stack_term("negative-zero", -0.0, 1.0),
                stack_term("other", 1.0, 1.0),
            ],
        ),
        Err(CorrelatedStackError::InvalidTermField {
            index: 0,
            field: "signed_sensitivity",
            issue: ScalarIssue::NonCanonicalNegativeZero,
            ..
        })
    ));
    for deviation in [0.0, -1.0, f64::NAN] {
        assert!(matches!(
            propagate_correlated_stack(
                &model,
                &[
                    stack_term("bad-deviation", 1.0, deviation),
                    stack_term("other", 1.0, 1.0),
                ],
            ),
            Err(CorrelatedStackError::InvalidTermField {
                index: 0,
                field: "standard_deviation",
                ..
            })
        ));
    }
    assert_eq!(
        propagate_correlated_stack(
            &model,
            &[
                stack_term("overflow", f64::MAX, 2.0),
                stack_term("other", 1.0, 1.0),
            ],
        ),
        Err(CorrelatedStackError::InvalidDerived {
            quantity: CorrelatedDerivedQuantity::ScaledSensitivity,
            term_index: Some(0),
            issue: ScalarIssue::NonFinite,
        })
    );
    assert_eq!(
        propagate_correlated_stack(
            &model,
            &[
                stack_term("underflow", f64::MIN_POSITIVE, f64::MIN_POSITIVE),
                stack_term("other", 1.0, 1.0),
            ],
        ),
        Err(CorrelatedStackError::InvalidDerived {
            quantity: CorrelatedDerivedQuantity::ScaledSensitivity,
            term_index: Some(0),
            issue: ScalarIssue::Underflow,
        })
    );
}

#[test]
fn correlated_stack_refuses_a_false_zero_from_normalization_underflow() {
    let model = AdmittedCorrelationModel::try_new(
        "gear/process-runout",
        NonZeroU64::new(1).expect("one is nonzero"),
        [0x5a; 32],
        3,
        vec![1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0],
    )
    .expect("singular manufactured factor is structurally admissible");
    let terms = [
        stack_term("first", 4.0, 1.0),
        stack_term("opposite", -4.0, 1.0),
        stack_term("subnormal-residual", f64::from_bits(1), 1.0),
    ];

    assert_eq!(
        propagate_correlated_stack(&model, &terms),
        Err(CorrelatedStackError::InvalidDerived {
            quantity: CorrelatedDerivedQuantity::NormalizedSensitivity,
            term_index: Some(2),
            issue: ScalarIssue::Underflow,
        })
    );
}

#[test]
fn exact_zero_sensitivities_publish_zero_only_with_bound_model_and_axes() {
    let model = correlation_model(0.8, 0.6);
    let terms = [
        stack_term("carrier-runout", 0.0, 1.0),
        stack_term("gear-eccentricity", 0.0, 2.0),
    ];
    let receipt = propagate_correlated_stack(&model, &terms).expect("exact zeros are explicit");
    assert_eq!(receipt.independent_standard_deviation().to_bits(), 0);
    assert_eq!(receipt.independent_variance().to_bits(), 0);
    assert_eq!(receipt.correlated_standard_deviation().to_bits(), 0);
    assert_eq!(receipt.correlated_variance().to_bits(), 0);
    assert_eq!(receipt.correlation_variance_delta().to_bits(), 0);
    assert_eq!(receipt.model(), &model);
    assert_eq!(receipt.terms(), terms);
}

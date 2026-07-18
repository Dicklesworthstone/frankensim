//! G0/G3/G5 conformance for dependency-aware grouped tolerance allocation.

use std::num::NonZeroU64;

use fs_toleralloc::{
    Action, ColorRank, DependencyGroup, DependencyGroupId, Feature, GroupedAllocationError,
    GroupedAllocationResource, GroupedDependenceIdentity, GroupedDependenceModel,
    GroupedDerivedQuantity, GroupedFeature, ScalarIssue, allocate, allocate_grouped,
};

fn nz(value: u64) -> NonZeroU64 {
    NonZeroU64::new(value).expect("fixture versions are nonzero")
}

fn feature(name: &str, sensitivity: f64, cost_coeff: f64) -> Feature {
    Feature {
        name: name.into(),
        sensitivity,
        sensitivity_color: ColorRank::Verified,
        cost_coeff,
        baseline_tolerance: 1.0,
    }
}

fn grouped_model() -> GroupedDependenceModel {
    GroupedDependenceModel {
        identity: GroupedDependenceIdentity::try_new(
            "gear/coherent-tolerance-groups",
            nz(1),
            [0x73; 32],
        )
        .expect("fixture identity is canonical"),
        groups: vec![
            DependencyGroup {
                key: "shaft-stack".into(),
            },
            DependencyGroup {
                key: "bearing-seat".into(),
            },
        ],
        features: vec![
            GroupedFeature {
                feature: feature("shaft-a", 1.0, 1.0),
                group: DependencyGroupId(0),
            },
            GroupedFeature {
                feature: feature("shaft-b", 1.0, 1.0),
                group: DependencyGroupId(0),
            },
            GroupedFeature {
                feature: feature("bearing-seat", 2.0, 2.0),
                group: DependencyGroupId(1),
            },
        ],
    }
}

fn assert_close(actual: f64, expected: f64) {
    let scale = actual.abs().max(expected.abs()).max(1.0);
    assert!(
        (actual - expected).abs() <= 128.0 * f64::EPSILON * scale,
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}

#[test]
fn perfect_group_allocation_matches_the_closed_form_oracle() {
    let model = grouped_model();
    let receipt = allocate_grouped(&model, 2.0, 1.0).expect("grouped allocation evaluates");

    assert_eq!(receipt.schema_version(), 1);
    assert_eq!(receipt.model(), &model);
    assert_eq!(receipt.variance_budget(), 2.0);
    assert_eq!(receipt.tolerance_to_sigma(), 1.0);
    assert_eq!(receipt.items().len(), 3);
    for (index, item) in receipt.items().iter().enumerate() {
        assert_eq!(item.feature_index(), index);
        assert_close(item.tolerance(), 0.5);
        assert_eq!(item.action(), Action::Tighten);
    }
    assert_eq!(receipt.items()[0].group(), DependencyGroupId(0));
    assert_eq!(receipt.items()[1].group(), DependencyGroupId(0));
    assert_eq!(receipt.items()[2].group(), DependencyGroupId(1));
    assert_close(receipt.items()[0].cost_contribution(), 2.0);
    assert_close(receipt.items()[1].cost_contribution(), 2.0);
    assert_close(receipt.items()[2].cost_contribution(), 4.0);
    assert_close(receipt.items()[0].standard_deviation_loading(), 0.5);
    assert_close(receipt.items()[1].standard_deviation_loading(), 0.5);
    assert_close(receipt.items()[2].standard_deviation_loading(), 1.0);

    assert_eq!(receipt.groups().len(), 2);
    let shaft = &receipt.groups()[0];
    assert_eq!(shaft.group(), DependencyGroupId(0));
    assert_eq!(shaft.feature_count(), 2);
    assert_close(shaft.standard_deviation(), 1.0);
    assert_close(shaft.variance(), 1.0);
    assert_close(shaft.independent_variance(), 0.5);
    assert_close(shaft.dependency_variance_delta(), 0.5);
    assert_close(shaft.total_cost(), 4.0);

    let bearing = &receipt.groups()[1];
    assert_eq!(bearing.group(), DependencyGroupId(1));
    assert_eq!(bearing.feature_count(), 1);
    assert_close(bearing.standard_deviation(), 1.0);
    assert_close(bearing.variance(), 1.0);
    assert_close(bearing.independent_variance(), 1.0);
    assert_close(bearing.dependency_variance_delta(), 0.0);
    assert_close(bearing.total_cost(), 4.0);

    // Independent arithmetic oracle from the published tolerances; this does
    // not trust the receipt's achieved-variance field.
    let tolerances = receipt
        .items()
        .iter()
        .map(|item| item.tolerance())
        .collect::<Vec<_>>();
    let shaft_loading = tolerances[0] + tolerances[1];
    let bearing_loading = 2.0 * tolerances[2];
    let grouped_variance = shaft_loading * shaft_loading + bearing_loading * bearing_loading;
    assert_close(grouped_variance, 2.0);
    assert_close(receipt.achieved_variance(), 2.0);
    assert_close(receipt.budget_residual(), 0.0);
    assert_close(receipt.independent_variance(), 1.5);
    assert_close(receipt.dependency_variance_delta(), 0.5);
    assert_close(receipt.total_cost(), 8.0);
    assert_close(receipt.closed_form_cost(), 8.0);
    assert_close(receipt.cost_residual(), 0.0);
    assert_close(receipt.log_scale_correction(), 0.0);
    assert_close(receipt.max_stationarity_log_residual(), 0.0);
}

#[test]
fn unequal_features_in_one_dependency_group_match_the_closed_form_oracle() {
    let model = GroupedDependenceModel {
        identity: GroupedDependenceIdentity::try_new(
            "gear/unequal-coherent-pair",
            nz(1),
            [0x41; 32],
        )
        .expect("fixture identity is canonical"),
        groups: vec![DependencyGroup {
            key: "coherent-pair".into(),
        }],
        features: vec![
            GroupedFeature {
                feature: feature("first", 1.0, 1.0),
                group: DependencyGroupId(0),
            },
            GroupedFeature {
                feature: feature("second", 4.0, 1.0),
                group: DependencyGroupId(0),
            },
        ],
    };

    let receipt = allocate_grouped(&model, 1.0, 1.0).expect("grouped allocation evaluates");

    assert_close(receipt.items()[0].tolerance(), 1.0 / 3.0);
    assert_close(receipt.items()[1].tolerance(), 1.0 / 6.0);
    assert_close(receipt.items()[0].standard_deviation_loading(), 1.0 / 3.0);
    assert_close(receipt.items()[1].standard_deviation_loading(), 2.0 / 3.0);
    assert_close(receipt.items()[0].cost_contribution(), 3.0);
    assert_close(receipt.items()[1].cost_contribution(), 6.0);
    assert_eq!(receipt.items()[0].action(), Action::Tighten);
    assert_eq!(receipt.items()[1].action(), Action::Tighten);

    assert_eq!(receipt.groups().len(), 1);
    let published_group_loading = receipt.items()[0].standard_deviation_loading()
        + receipt.items()[1].standard_deviation_loading();
    assert_close(
        receipt.groups()[0].standard_deviation(),
        published_group_loading,
    );
    assert_close(receipt.groups()[0].standard_deviation(), 1.0);
    assert_close(receipt.groups()[0].variance(), 1.0);
    assert_close(receipt.achieved_variance(), 1.0);
    assert_close(receipt.budget_residual(), 0.0);
    assert_close(receipt.total_cost(), 9.0);
    assert_close(receipt.closed_form_cost(), 9.0);
    assert_close(receipt.max_stationarity_log_residual(), 0.0);
}

#[test]
fn legacy_independent_allocation_violates_the_true_group_budget() {
    let model = grouped_model();
    let plain_features = model
        .features
        .iter()
        .map(|grouped| grouped.feature.clone())
        .collect::<Vec<_>>();
    let legacy = allocate(&plain_features, 2.0, 1.0).expect("legacy allocation evaluates");
    assert_close(legacy.achieved_variance, 2.0);

    let shaft_loading = legacy.items[0].tolerance + legacy.items[1].tolerance;
    let bearing_loading = 2.0 * legacy.items[2].tolerance;
    let true_grouped_variance = shaft_loading * shaft_loading + bearing_loading * bearing_loading;
    assert!(
        true_grouped_variance > 2.8,
        "independence should materially violate the grouped budget: {true_grouped_variance}"
    );

    let grouped = allocate_grouped(&model, 2.0, 1.0).expect("grouped allocation evaluates");
    let feasible_scale = (2.0 / true_grouped_variance).sqrt();
    let feasible_legacy_cost = legacy
        .items
        .iter()
        .zip(&plain_features)
        .map(|(item, feature)| feature.cost_coeff / (item.tolerance * feasible_scale))
        .sum::<f64>();
    assert!(
        grouped.total_cost() < feasible_legacy_cost,
        "closed-form grouped optimum must beat a normalized legacy candidate"
    );
}

#[test]
fn g3_common_sensitivity_and_cost_rescaling_is_covariant() {
    let model = grouped_model();
    let base = allocate_grouped(&model, 2.0, 1.0).expect("base evaluates");

    let mut sensitivity_scaled = model.clone();
    for grouped in &mut sensitivity_scaled.features {
        grouped.feature.sensitivity *= 4.0;
    }
    let scaled = allocate_grouped(&sensitivity_scaled, 32.0, 1.0).expect("scaled model evaluates");
    for (left, right) in base.items().iter().zip(scaled.items()) {
        assert_close(left.tolerance(), right.tolerance());
        assert_eq!(left.action(), right.action());
    }
    assert_close(scaled.total_cost(), base.total_cost());

    let mut cost_scaled = model.clone();
    for grouped in &mut cost_scaled.features {
        grouped.feature.cost_coeff *= 9.0;
    }
    let scaled = allocate_grouped(&cost_scaled, 2.0, 1.0).expect("cost scaling evaluates");
    for (left, right) in base.items().iter().zip(scaled.items()) {
        assert_close(left.tolerance(), right.tolerance());
        assert_eq!(left.action(), right.action());
    }
    assert_close(scaled.total_cost(), 9.0 * base.total_cost());
    assert_close(scaled.achieved_variance(), base.achieved_variance());

    let mut renamed = model;
    renamed.groups[0].key = "shaft-stack-renamed".into();
    let renamed_receipt =
        allocate_grouped(&renamed, 2.0, 1.0).expect("renamed provenance evaluates");
    assert_close(renamed_receipt.total_cost(), base.total_cost());
    assert_close(
        renamed_receipt.achieved_variance(),
        base.achieved_variance(),
    );
    assert_ne!(renamed_receipt.model(), base.model());
}

#[test]
fn singleton_groups_reduce_to_the_independent_allocator() {
    let identity = GroupedDependenceIdentity::try_new("gear/singleton-groups", nz(1), [0x29; 32])
        .expect("identity is canonical");
    let features = vec![feature("first", 1.5, 2.0), feature("second", 0.75, 3.0)];
    let model = GroupedDependenceModel {
        identity,
        groups: vec![
            DependencyGroup {
                key: "first-group".into(),
            },
            DependencyGroup {
                key: "second-group".into(),
            },
        ],
        features: features
            .iter()
            .cloned()
            .enumerate()
            .map(|(index, feature)| GroupedFeature {
                feature,
                group: DependencyGroupId(index as u16),
            })
            .collect(),
    };
    let independent = allocate(&features, 0.75, 3.0).expect("independent allocation evaluates");
    let grouped = allocate_grouped(&model, 0.75, 3.0).expect("singleton groups evaluate");
    for (plain, dependent) in independent.items.iter().zip(grouped.items()) {
        assert_close(plain.tolerance, dependent.tolerance());
        assert_eq!(plain.action, dependent.action());
    }
    assert_close(independent.total_cost, grouped.total_cost());
    assert_close(independent.achieved_variance, grouped.achieved_variance());
    assert_close(grouped.dependency_variance_delta(), 0.0);
}

#[test]
fn g5_identical_input_replays_the_complete_receipt() {
    let model = grouped_model();
    let first = allocate_grouped(&model, 2.0, 1.0).expect("first evaluates");
    let second = allocate_grouped(&model, 2.0, 1.0).expect("second evaluates");
    assert_eq!(first, second);
    assert_eq!(
        first.model().identity.namespace(),
        "gear/coherent-tolerance-groups"
    );
    assert_eq!(first.model().identity.schema_version().get(), 1);
    assert_eq!(first.model().identity.semantic_digest(), [0x73; 32]);
}

#[test]
fn grouped_model_admission_refuses_ambiguous_or_incomplete_dependencies() {
    let mut model = grouped_model();
    model.groups.clear();
    assert_eq!(
        allocate_grouped(&model, 2.0, 1.0),
        Err(GroupedAllocationError::ResourceLimit {
            resource: GroupedAllocationResource::Groups,
            actual: 0,
            max: 128,
        })
    );

    let mut model = grouped_model();
    model.features.clear();
    assert_eq!(
        allocate_grouped(&model, 2.0, 1.0),
        Err(GroupedAllocationError::ResourceLimit {
            resource: GroupedAllocationResource::Features,
            actual: 0,
            max: 128,
        })
    );

    let mut model = grouped_model();
    model.groups[0].key = "Shaft-Stack".into();
    assert!(matches!(
        allocate_grouped(&model, 2.0, 1.0),
        Err(GroupedAllocationError::InvalidGroupKey { index: 0, .. })
    ));

    let mut model = grouped_model();
    model.groups[1].key = "shaft-stack".into();
    assert!(matches!(
        allocate_grouped(&model, 2.0, 1.0),
        Err(GroupedAllocationError::DuplicateGroupKey {
            first_index: 0,
            duplicate_index: 1,
            ..
        })
    ));

    let mut model = grouped_model();
    model.features[1].feature.name = "SHAFT-A".into();
    assert!(matches!(
        allocate_grouped(&model, 2.0, 1.0),
        Err(GroupedAllocationError::DuplicateFeatureName {
            first_index: 0,
            duplicate_index: 1,
            ..
        })
    ));

    let mut model = grouped_model();
    model.features[0].group = DependencyGroupId(7);
    assert_eq!(
        allocate_grouped(&model, 2.0, 1.0),
        Err(GroupedAllocationError::InvalidGroupReference {
            feature_index: 0,
            group: DependencyGroupId(7),
            available: 2,
        })
    );

    let mut model = grouped_model();
    model
        .features
        .retain(|feature| feature.group == DependencyGroupId(0));
    assert_eq!(
        allocate_grouped(&model, 2.0, 1.0),
        Err(GroupedAllocationError::EmptyGroup {
            group: DependencyGroupId(1),
            key: "bearing-seat".into(),
        })
    );
}

#[test]
fn grouped_allocation_refuses_bad_scalars_and_lost_log_contributions() {
    for (budget, k, argument) in [
        (0.0, 1.0, "variance_budget"),
        (f64::NAN, 1.0, "variance_budget"),
        (2.0, 0.0, "k"),
        (2.0, f64::INFINITY, "k"),
    ] {
        assert!(matches!(
            allocate_grouped(&grouped_model(), budget, k),
            Err(GroupedAllocationError::InvalidArgument {
                argument: actual,
                ..
            }) if actual == argument
        ));
    }

    let mut model = grouped_model();
    model.features[0].feature.sensitivity = -0.0;
    assert_eq!(
        allocate_grouped(&model, 2.0, 1.0),
        Err(GroupedAllocationError::InvalidFeatureField {
            index: 0,
            feature: "shaft-a".into(),
            field: "sensitivity",
            issue: ScalarIssue::NonPositive,
        })
    );

    let mut model = grouped_model();
    model.features[0].feature.sensitivity = f64::MAX;
    model.features[0].feature.cost_coeff = f64::MAX;
    model.features[1].feature.sensitivity = f64::MIN_POSITIVE;
    model.features[1].feature.cost_coeff = f64::MIN_POSITIVE;
    assert!(matches!(
        allocate_grouped(&model, 2.0, 1.0),
        Err(GroupedAllocationError::InvalidDerived {
            quantity: GroupedDerivedQuantity::LogSumExpContribution,
            group: Some(DependencyGroupId(0)),
            issue: ScalarIssue::Underflow,
            ..
        })
    ));
}

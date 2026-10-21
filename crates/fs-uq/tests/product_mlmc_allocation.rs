use fs_uq::{allocate_mlmc, MlmcAllocationError, MlmcAllocationPolicy};

fn policy(target: f64, calls: usize) -> MlmcAllocationPolicy {
    MlmcAllocationPolicy { target_standard_error: target, minimum_samples_per_level: 2,
        correction_std_dev_floor: 1e-12, maximum_model_calls: calls }
}

#[test]
fn analytic_giles_allocation_prices_both_halves_and_meets_predicted_target() {
    let result = allocate_mlmc(&[4.0, 2.0, 1.0], &[1.0, 3.0, 12.0], policy(1.0, 1000)).unwrap();
    assert_eq!(result.samples_per_level, [48, 12, 4]);
    assert_eq!(result.model_calls, 80);
    assert_eq!(result.declared_cost, Some(156.0));
    let expected = fs_math::det::sqrt(16.0 / 48.0 + 4.0 / 12.0 + 1.0 / 4.0);
    assert!((result.predicted_standard_error - expected).abs() < 1e-14);
    assert!(result.predicted_standard_error <= 1.0);
    assert!(result.samples_per_level[0] > 10 * result.samples_per_level[2]);
}

#[test]
fn units_and_common_cost_scale_do_not_change_the_allocation() {
    let reference = allocate_mlmc(&[4.0, 2.0, 1.0], &[1.0, 3.0, 12.0], policy(1.0, 1000)).unwrap();
    for scale in [1e-150, 1e150] {
        let mut p = policy(scale, 1000);
        p.correction_std_dev_floor *= scale;
        let scaled = allocate_mlmc(&[4.0 * scale, 2.0 * scale, scale], &[1.0, 3.0, 12.0], p).unwrap();
        assert_eq!(scaled.samples_per_level, reference.samples_per_level);
        assert!((scaled.predicted_standard_error / scale - reference.predicted_standard_error).abs() < 1e-13);
        let priced = allocate_mlmc(&[4.0, 2.0, 1.0], &[scale, 3.0 * scale, 12.0 * scale], policy(1.0, 1000)).unwrap();
        assert_eq!(priced.samples_per_level, reference.samples_per_level);
        assert!((priced.declared_cost.unwrap() / scale - 156.0).abs() < 1e-10);
    }
}

#[test]
fn constant_pilot_uses_explicit_floor_and_still_allocates_every_level() {
    let mut p = policy(0.1, 1000);
    p.correction_std_dev_floor = 0.25;
    let result = allocate_mlmc(&[0.0, 0.0], &[1.0, 1.0], p).unwrap();
    assert_eq!(result.allocation_standard_deviations, [0.25, 0.25]);
    assert!(result.samples_per_level.iter().all(|&n| n > 2));
    assert!(result.predicted_standard_error > 0.0 && result.predicted_standard_error <= 0.1);
    let tiny = allocate_mlmc(&[0.0; 3], &[1.0; 3], policy(1.0, 10)).unwrap();
    assert_eq!(tiny.samples_per_level, [2, 2, 2]);
    assert_eq!(tiny.model_calls, 10);
    assert!(tiny.predicted_standard_error > 0.0);
}

#[test]
fn unavailable_resources_or_malformed_data_never_publish_an_underallocation() {
    assert!(matches!(allocate_mlmc(&[4.0, 2.0, 1.0], &[1.0, 3.0, 12.0], policy(1.0, 79)),
        Err(MlmcAllocationError::BudgetExceeded { .. })));
    for p in [policy(0.0, 100), policy(f64::NAN, 100), policy(1.0, 3), policy(1.0, 1_000_001),
        MlmcAllocationPolicy { minimum_samples_per_level: 1, ..policy(1.0, 100) },
        MlmcAllocationPolicy { correction_std_dev_floor: 0.0, ..policy(1.0, 100) }] {
        assert!(allocate_mlmc(&[1.0, 1.0], &[1.0, 2.0], p).is_err());
    }
    for (deviations, costs) in [(vec![], vec![]), (vec![1.0], vec![1.0, 2.0]),
        (vec![f64::NAN], vec![1.0]), (vec![-1.0], vec![1.0]),
        (vec![1.0], vec![0.0]), (vec![1.0], vec![f64::INFINITY]),
        (vec![1.0; 65], vec![1.0; 65])] {
        assert!(allocate_mlmc(&deviations, &costs, policy(1.0, 1000)).is_err());
    }
    assert!(allocate_mlmc(&[f64::MAX], &[1.0], policy(f64::MIN_POSITIVE, 1000)).is_err());
}

#[test]
fn dimensional_cost_overflow_does_not_invalidate_finite_counts_or_error() {
    let result = allocate_mlmc(&[1.0], &[f64::MAX], policy(1.0, 4)).unwrap();
    assert_eq!(result.samples_per_level, [2]);
    assert_eq!(result.model_calls, 2);
    assert_eq!(result.declared_cost, None);
    assert!(result.predicted_standard_error.is_finite());
    let mut p = policy(1e-300, 1000);
    p.correction_std_dev_floor = 1e-310;
    let result = allocate_mlmc(&[2e-300], &[1.0], p).unwrap();
    assert!(result.predicted_standard_error > 0.0 && result.predicted_standard_error <= 1e-300);
}

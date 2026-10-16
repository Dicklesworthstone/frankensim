use fs_couple::iqn_ils::{IqnIls, IqnIlsConfig, IqnIlsError, MAX_HISTORY};

fn accelerator(dimension: usize) -> IqnIls {
    IqnIls::new(dimension, IqnIlsConfig::default()).unwrap()
}

fn assert_close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (&actual, &expected) in actual.iter().zip(expected) {
        assert!((actual - expected).abs() < 1.0e-10, "{actual} != {expected}");
    }
}

#[test]
fn opposing_residuals_do_not_disappear_into_a_scalar_mean() {
    let mut iqn = accelerator(2);
    let mut x: Vec<f64> = vec![0.0, 0.0];
    let mut plain = x.clone();
    for _ in 0..5 {
        // The signed mean is zero at every nonconverged iterate. Plain
        // staggering diverges; the exact fixed point is [1, -1].
        let image = vec![3.0 - 2.0 * x[0], -3.0 - 2.0 * x[1]];
        assert!(((image[0] - x[0]) + (image[1] - x[1])).abs() < 1.0e-12);
        x = iqn.step(&x, &image, 0.5).unwrap().values;
        plain = vec![3.0 - 2.0 * plain[0], -3.0 - 2.0 * plain[1]];
    }
    assert_close(&x, &[1.0, -1.0]);
    assert!((plain[0] - 1.0).abs() > 10.0);
}

#[test]
fn independent_modes_recover_the_coupled_linear_solution() {
    let mut iqn = accelerator(2);
    let mut x: Vec<f64> = vec![0.0, 0.0];
    let mut saw_full_rank = false;
    for _ in 0..8 {
        // (I-A)x=b gives x=[1,-2], independently by direct substitution.
        let image = vec![5.0 - 2.0 * x[0] + x[1], -1.5 + 0.5 * x[0] + 0.5 * x[1]];
        let step = iqn.step(&x, &image, 0.25).unwrap();
        saw_full_rank |= step.used_columns == 2;
        x = step.values;
    }
    assert!(saw_full_rank);
    assert_close(&x, &[1.0, -2.0]);
}

#[test]
fn startup_and_repeated_residual_use_the_declared_fallback() {
    let mut iqn = accelerator(2);
    let first = iqn.step(&[0.0, 0.0], &[2.0, -4.0], 0.25).unwrap();
    assert_eq!(first.values, vec![0.5, -1.0]);
    assert_eq!(first.used_columns, 0);
    assert_eq!(first.relaxation_omega, 0.25);
    let second = iqn.step(&first.values, &[2.5, -5.0], 0.25).unwrap();
    assert_eq!(second.values, vec![1.0, -2.0]);
    assert_eq!(second.used_columns, 0);
    assert_eq!(iqn.history_len(), 0);
}

#[test]
fn dependent_history_is_filtered_without_a_singular_solve() {
    let mut iqn = accelerator(2);
    iqn.step(&[0.0, 0.0], &[1.0, -1.0], 0.5).unwrap();
    for k in 1..20 {
        let k = f64::from(k);
        let step = iqn
            .step(&[k, -k], &[2.0 * k + 1.0, -2.0 * k - 1.0], 0.5)
            .unwrap();
        assert_eq!(step.used_columns, 1);
        assert!(step.values.iter().all(|value| value.is_finite()));
        assert!(iqn.history_len() <= IqnIlsConfig::default().max_history);
    }
}

#[test]
fn nearly_dependent_history_respects_the_rank_threshold() {
    let mut iqn = IqnIls::new(
        2,
        IqnIlsConfig {
            max_history: 4,
            relative_rank_tolerance: 1.0e-6,
        },
    )
    .unwrap();
    iqn.step(&[0.0, 0.0], &[0.0, 0.0], 0.5).unwrap();
    iqn.step(&[0.0, 0.0], &[1.0, 1.0], 0.5).unwrap();
    let step = iqn
        .step(&[0.0, 0.0], &[2.0, 2.0 + 1.0e-10], 0.5)
        .unwrap();
    assert_eq!(step.used_columns, 1);
    assert!(step.values.iter().all(|value| value.is_finite()));
}

#[test]
fn malformed_input_is_transactional_and_retryable() {
    let mut iqn = accelerator(2);
    iqn.step(&[0.0, 0.0], &[1.0, 2.0], 0.5).unwrap();
    let snapshot = iqn.clone();
    assert!(matches!(
        iqn.step(&[0.0], &[1.0, 2.0], 0.5),
        Err(IqnIlsError::DimensionMismatch { .. })
    ));
    assert!(iqn.step(&[0.0, 0.0], &[f64::NAN, 1.0], 0.5).is_err());
    assert!(
        iqn.step(&[f64::MAX, 0.0], &[-f64::MAX, 1.0], 0.5)
            .is_err()
    );
    assert!(
        iqn.step(&[0.0, 0.0], &[1.0, 2.0], f64::INFINITY)
            .is_err()
    );
    assert_eq!(iqn, snapshot);
    let mut clean = snapshot;
    assert_eq!(
        iqn.step(&[0.5, 1.0], &[0.0, 1.5], 0.5),
        clean.step(&[0.5, 1.0], &[0.0, 1.5], 0.5)
    );
}

#[test]
fn output_overflow_does_not_commit_a_sample() {
    let mut iqn = accelerator(1);
    let snapshot = iqn.clone();
    assert!(iqn.step(&[0.0], &[f64::MAX], 2.0).is_err());
    assert_eq!(iqn, snapshot);
    assert_eq!(iqn.step(&[1.0], &[1.0], 0.5).unwrap().values, vec![1.0]);
}

#[test]
fn full_state_checkpoint_replays_and_reset_removes_old_map_history() {
    let mut iqn = accelerator(2);
    iqn.step(&[0.0, 0.0], &[1.0, 2.0], 0.5).unwrap();
    iqn.step(&[0.5, 1.0], &[0.0, 1.5], 0.5).unwrap();
    let mut restored = iqn.clone();
    assert_eq!(
        iqn.step(&[0.3, 1.3], &[0.4, 1.35], 0.5),
        restored.step(&[0.3, 1.3], &[0.4, 1.35], 0.5)
    );
    iqn.reset();
    assert_eq!(iqn, accelerator(2));
}

#[test]
fn configuration_refuses_invalid_dimensions_and_budgets() {
    assert!(IqnIls::new(0, IqnIlsConfig::default()).is_err());
    for max_history in [0, MAX_HISTORY + 1, usize::MAX] {
        assert!(
            IqnIls::new(
                2,
                IqnIlsConfig {
                    max_history,
                    ..IqnIlsConfig::default()
                }
            )
            .is_err()
        );
    }
    for relative_rank_tolerance in [0.0, -1.0, 1.0, f64::NAN, f64::INFINITY] {
        assert!(
            IqnIls::new(
                2,
                IqnIlsConfig {
                    relative_rank_tolerance,
                    ..IqnIlsConfig::default()
                }
            )
            .is_err()
        );
    }
}

#[test]
fn fallback_preserves_signed_and_zero_aitken_factors() {
    for omega in [-3.0, -0.5, 0.0, 3.0] {
        let mut iqn = accelerator(2);
        let step = iqn.step(&[2.0, -1.0], &[3.0, 1.0], omega).unwrap();
        assert_eq!(step.values, vec![2.0 + omega, -1.0 + 2.0 * omega]);
        assert_eq!(step.used_columns, 0);
        assert_eq!(step.relaxation_omega, omega);
    }
}

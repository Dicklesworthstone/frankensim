//! G0/G3/G4 conditional moments for independently noisy threshold sensors.

use fs_assimilate::{
    AssimError, BatchOutcome, Belief, CensorDirection, Observation, ObservationBatch,
    ObservationDisposition, RobustAssimilation, RobustObservation, assimilate_all,
    assimilate_observation_batch, point_sensor,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};

fn with_cx<R>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 713,
                kernel_id: 91,
                tile: 0,
                iteration: 0,
            },
            budget,
            ExecMode::Deterministic,
        )
        .with_time_source(&clock))
    });
    assert!(pool.stats().quiescent());
    result
}

fn censor(threshold: f64, noise: f64, direction: CensorDirection, id: &str) -> RobustObservation {
    RobustObservation::censored(
        point_sensor(0, 1, 123.0, noise, id).unwrap(),
        threshold,
        direction,
    )
    .unwrap()
}

fn update(prior: &Belief, records: Vec<RobustObservation>, cx: &Cx<'_>) -> RobustAssimilation {
    let batch = ObservationBatch::independent(records, cx).unwrap();
    assimilate_observation_batch(prior, &batch, cx).unwrap()
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0e-300),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn a_threshold_is_an_event_not_a_point_reading_in_both_directions() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let prior = Belief::scalar(0.0, 1.0).unwrap();
        for (direction, sign) in [
            (CensorDirection::AtOrAbove, 1.0),
            (CensorDirection::AtOrBelow, -1.0),
        ] {
            let result = update(&prior, vec![censor(0.0, 1.0, direction, "threshold")], cx);
            let posterior = result.posterior().unwrap();
            close(
                posterior.mean()[0],
                sign / core::f64::consts::PI.sqrt(),
                4.0e-15,
            );
            close(
                posterior.variance(0).unwrap(),
                1.0 - 1.0 / core::f64::consts::PI,
                4.0e-15,
            );
            assert_eq!(result.audit().outcome(), BatchOutcome::Updated);
            assert_eq!(result.audit().effective_dof(), 1);
            assert_eq!(
                result.audit().entries()[0].1,
                ObservationDisposition::AcceptedCensoredMoments
            );
            assert!(!result.audit().uses_sequential_moment_approximation());
            assert!(
                result
                    .audit()
                    .receipt_id()
                    .starts_with("robust-observation-audit:v2:")
            );
        }
    });
}

#[test]
fn conditional_cross_covariance_matches_an_independent_gaussian_integral() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let prior = Belief::new(vec![0.3, -0.2], vec![vec![2.0, 0.5], vec![0.5, 1.0]], cx).unwrap();
        let reading = Observation::new(vec![1.0, 0.5], -900.0, 0.25, "linear-combination").unwrap();
        let event = RobustObservation::censored(reading, 1.2, CensorDirection::AtOrAbove).unwrap();
        let result = update(&prior, vec![event], cx);
        let posterior = result.posterior().unwrap();
        // Independent 80-digit erfc evaluation of the scalar normal tail,
        // followed by Gaussian conditioning (not this implementation's CF).
        let expected_mean = [1.856_428_945_886_533_5, 0.491_746_198_171_792_66];
        let expected_cov = [
            [0.744_850_645_821_434_3, -0.057_844_157_412_695_86],
            [-0.057_844_157_412_695_86, 0.752_069_263_372_135_2],
        ];
        for row in 0..2 {
            close(posterior.mean()[row], expected_mean[row], 2.0e-14);
            assert!(posterior.covariance()[row][row] < prior.covariance()[row][row]);
            for column in 0..2 {
                close(
                    posterior.covariance()[row][column],
                    expected_cov[row][column],
                    8.0e-14,
                );
            }
        }
        posterior.validate(cx).unwrap();
    });
}

#[test]
fn normal_tail_oracles_cover_central_switch_and_underflowing_tail_probability() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        // Y has unit predictive variance; X has variance 0.8 and noise 0.2.
        // Reference tail moments come from a 90-digit erfc evaluation.
        let cases = [
            (-4.0, 0.000_133_834_464_468_575_14, 0.999_464_644_230_461_8),
            (-1.0, 0.287_599_970_939_178_36, 0.629_686_285_776_605_4),
            (1.499_999, 1.938_676_316_169_178_3, 0.149_546_676_702_155_46),
            (1.5, 1.938_677_166_622_543_2, 0.149_546_593_550_202_7),
            (8.0, 8.121_368_112_236_113, 0.014_324_883_443_340_91),
            (40.0, 40.024_968_847_207_264, 0.000_622_668_378_591_388_8),
            (100.0, 100.009_998_000_999_26, 0.000_099_940_049_948_263_45),
        ];
        let prior = Belief::scalar(0.0, 0.8).unwrap();
        for (threshold, mean, variance) in cases {
            for (direction, sign) in [
                (CensorDirection::AtOrAbove, 1.0),
                (CensorDirection::AtOrBelow, -1.0),
            ] {
                let result = update(
                    &prior,
                    vec![censor(sign * threshold, 0.2, direction, "tail")],
                    cx,
                );
                let posterior = result.posterior().unwrap();
                close(posterior.mean()[0], sign * 0.8 * mean, 8.0e-14);
                close(
                    posterior.covariance()[0][0],
                    0.16 + 0.64 * variance,
                    8.0e-14,
                );
            }
        }
        // A direct 1 - lambda*(lambda-a) loses this variance completely.
        let prior = Belief::scalar(0.0, 1.0).unwrap();
        let result = update(
            &prior,
            vec![censor(1.0e8, 1.0e-20, CensorDirection::AtOrAbove, "rare")],
            cx,
        );
        let posterior = result.posterior().unwrap();
        close(posterior.mean()[0], 100_000_000.000_000_01, 3.0e-16);
        close(
            posterior.covariance()[0][0],
            9.999_999_999_999_995e-17 + 1.0e-20,
            2.0e-14,
        );
    });
}

#[test]
fn physical_covariance_survives_underflow_of_dimensionless_tail_variance() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let prior = Belief::scalar(0.0, 2.0_f64.powi(864)).unwrap();
        for (direction, sign) in [
            (CensorDirection::AtOrAbove, 1.0),
            (CensorDirection::AtOrBelow, -1.0),
        ] {
            let observation = Observation::new(
                vec![2.0_f64.powi(-266)],
                0.0,
                2.0_f64.powi(-996),
                "extreme-units",
            )
            .unwrap();
            let event =
                RobustObservation::censored(observation, sign * 2.0_f64.powi(730), direction)
                    .unwrap();
            let result = update(&prior, vec![event], cx);
            let posterior = result.posterior().unwrap();
            // a=2^564, so Var(Z|Z>=a)~2^-1128 is below the least subnormal.
            // Its standard deviation 2^-564 times Ph/sqrt(s)=2^432 remains
            // representable. The physical residual variance is 2^-264,
            // dwarfing the point-conditioned noise covariance 2^-464.
            close(posterior.mean()[0], sign * 2.0_f64.powi(996), 2.0e-15);
            close(posterior.covariance()[0][0], 2.0_f64.powi(-264), 2.0e-15);
            posterior.validate(cx).unwrap();
        }
    });
}

#[test]
fn clipping_equals_the_corresponding_tail_and_ignores_unobserved_placeholder_values() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let prior = Belief::scalar(0.0, 1.0).unwrap();
        for (endpoint, direction) in [
            (-3.0, CensorDirection::AtOrBelow),
            (5.0, CensorDirection::AtOrAbove),
        ] {
            let clipping = RobustObservation::saturated(
                point_sensor(0, 1, endpoint, 0.2, "clip").unwrap(),
                -3.0,
                5.0,
            )
            .unwrap();
            let clipped = update(&prior, vec![clipping], cx);
            let event = update(&prior, vec![censor(endpoint, 0.2, direction, "event")], cx);
            assert_eq!(clipped.posterior(), event.posterior());
            assert_eq!(
                clipped.audit().entries()[0].1,
                ObservationDisposition::AcceptedSaturatedMoments
            );
            let different_placeholder = RobustObservation::censored(
                point_sensor(0, 1, -2000.0, 0.2, "event").unwrap(),
                endpoint,
                direction,
            )
            .unwrap();
            let replay = update(&prior, vec![different_placeholder], cx);
            assert_eq!(event.posterior(), replay.posterior());
            assert_ne!(event.audit().receipt_id(), replay.audit().receipt_id());
        }
    });
}

#[test]
fn ordinary_readings_precede_events_and_multiple_events_are_explicit_approximations() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let prior = Belief::scalar(0.0, 1.0).unwrap();
        let ordinary = point_sensor(0, 1, 1.0, 1.0, "ordinary").unwrap();
        let event = censor(0.5, 0.5, CensorDirection::AtOrAbove, "a-event");
        let missing = RobustObservation::missing("missing").unwrap();
        let records = vec![
            event.clone(),
            missing.clone(),
            RobustObservation::available(ordinary.clone()),
        ];
        let single = update(&prior, records.clone(), cx);
        let posterior = single.posterior().unwrap();
        close(
            posterior.mean()[0],
            0.5 + 1.0 / (2.0 * core::f64::consts::PI).sqrt(),
            4.0e-15,
        );
        close(
            posterior.covariance()[0][0],
            0.5 - 0.5 / core::f64::consts::PI,
            4.0e-15,
        );
        assert_eq!(single.audit().effective_dof(), 2);
        assert!(!single.audit().uses_sequential_moment_approximation());
        let mut multiple = records;
        multiple.push(censor(1.5, 0.3, CensorDirection::AtOrBelow, "b-event"));
        let first = update(&prior, multiple.clone(), cx);
        multiple.reverse();
        let permuted = update(&prior, multiple, cx);
        assert_eq!(first.posterior(), permuted.posterior());
        assert!(first.audit().uses_sequential_moment_approximation());
        assert_eq!(first.audit().effective_dof(), 3);
        assert_eq!(
            first.audit().entries()[0].1,
            ObservationDisposition::AcceptedCensoredMomentApproximation
        );
        let just_ordinary = update(
            &prior,
            vec![missing, RobustObservation::available(ordinary.clone())],
            cx,
        );
        assert_eq!(
            just_ordinary.posterior(),
            Some(&assimilate_all(&prior, &[ordinary], cx).unwrap())
        );
    });
}

#[test]
fn measurement_unit_rescaling_and_almost_certain_events_preserve_moments() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let prior = Belief::scalar(2.0, 9.0).unwrap();
        let base = update(
            &prior,
            vec![censor(3.0, 4.0, CensorDirection::AtOrAbove, "base")],
            cx,
        );
        let scaled = RobustObservation::censored(
            Observation::new(vec![100.0], 0.0, 40_000.0, "scaled").unwrap(),
            300.0,
            CensorDirection::AtOrAbove,
        )
        .unwrap();
        let transformed = update(&prior, vec![scaled], cx);
        close(
            transformed.posterior().unwrap().mean()[0],
            base.posterior().unwrap().mean()[0],
            2.0e-14,
        );
        close(
            transformed.posterior().unwrap().covariance()[0][0],
            base.posterior().unwrap().covariance()[0][0],
            2.0e-14,
        );
        let certain = update(
            &prior,
            vec![censor(
                -100.0,
                4.0,
                CensorDirection::AtOrAbove,
                "almost-certain",
            )],
            cx,
        );
        close(
            certain.posterior().unwrap().mean()[0],
            prior.mean()[0],
            2.0e-15,
        );
        close(
            certain.posterior().unwrap().covariance()[0][0],
            prior.covariance()[0][0],
            2.0e-15,
        );
    });
}

#[test]
fn undeclared_dependence_delay_and_ambiguous_clipping_refuse_without_partial_updates() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let prior = Belief::scalar(0.0, 1.0).unwrap();
        let ordinary = point_sensor(0, 1, 1.0, 1.0, "ordinary").unwrap();
        let event = censor(0.0, 1.0, CensorDirection::AtOrAbove, "event");
        // A diagonal available-row covariance does not declare independence
        // of the censored noise from those rows or from other events.
        let legacy = ObservationBatch::new(
            vec![
                RobustObservation::available(ordinary.clone()),
                event.clone(),
            ],
            vec![vec![1.0]],
            cx,
        )
        .unwrap();
        let refused = assimilate_observation_batch(&prior, &legacy, cx).unwrap();
        assert_eq!(refused.audit().outcome(), BatchOutcome::RefusedPathology);
        assert!(refused.posterior().is_none());
        let correlated = ObservationBatch::new(
            vec![
                RobustObservation::available(ordinary.clone()),
                RobustObservation::available(point_sensor(0, 1, 0.5, 1.0, "correlated").unwrap()),
                event.clone(),
            ],
            vec![vec![1.0, 0.5], vec![0.5, 1.0]],
            cx,
        )
        .unwrap();
        let refused = assimilate_observation_batch(&prior, &correlated, cx).unwrap();
        assert_eq!(refused.audit().outcome(), BatchOutcome::RefusedPathology);
        assert!(refused.posterior().is_none());
        let delayed =
            RobustObservation::delayed(point_sensor(0, 1, 0.0, 1.0, "delayed").unwrap(), 0.5, 0.1)
                .unwrap();
        let ambiguous = RobustObservation::saturated(
            point_sensor(0, 1, 2.0, 1.0, "ambiguous").unwrap(),
            2.0,
            2.0,
        )
        .unwrap();
        for unsupported in [delayed, ambiguous] {
            let refused = update(
                &prior,
                vec![
                    RobustObservation::available(ordinary.clone()),
                    event.clone(),
                    unsupported,
                ],
                cx,
            );
            assert_eq!(refused.audit().outcome(), BatchOutcome::RefusedPathology);
            assert_eq!(refused.audit().effective_dof(), 0);
            assert!(refused.posterior().is_none());
        }
    });
}

#[test]
fn dimensional_mismatch_cancel_and_exhausted_budgets_never_publish_a_posterior() {
    let gate = CancelGate::new();
    let prior = Belief::scalar(0.0, 1.0).unwrap();
    let batch = with_cx(&gate, Budget::INFINITE, |cx| {
        let wrong = RobustObservation::censored(
            Observation::new(vec![1.0, 0.0], 0.0, 1.0, "wrong").unwrap(),
            0.0,
            CensorDirection::AtOrAbove,
        )
        .unwrap();
        let invalid = ObservationBatch::independent(vec![wrong], cx).unwrap();
        assert!(matches!(
            assimilate_observation_batch(&prior, &invalid, cx),
            Err(AssimError::DimMismatch { .. })
        ));
        ObservationBatch::independent(
            vec![censor(40.0, 1.0, CensorDirection::AtOrAbove, "bounded")],
            cx,
        )
        .unwrap()
    });
    with_cx(&gate, Budget::new().with_cost_quota(0), |cx| {
        assert!(matches!(
            assimilate_observation_batch(&prior, &batch, cx),
            Err(AssimError::BudgetRefused(_))
        ));
    });
    with_cx(
        &gate,
        Budget {
            poll_quota: 1,
            ..Budget::INFINITE
        },
        |cx| {
            assert!(matches!(
                assimilate_observation_batch(&prior, &batch, cx),
                Err(AssimError::BudgetRefused(_))
            ));
        },
    );
    gate.request();
    with_cx(&gate, Budget::INFINITE, |cx| {
        assert!(matches!(
            assimilate_observation_batch(&prior, &batch, cx),
            Err(AssimError::Cancelled { .. })
        ));
    });
}

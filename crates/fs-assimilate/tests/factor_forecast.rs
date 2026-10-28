//! Direct Gaussian matrix formulas are independent of the root-rotation kernel.
use fs_assimilate::factor::forecast::{LinearDynamics, NonlinearDynamics};
use fs_assimilate::factor::{FactorBelief, assimilate_scalar};
use fs_assimilate::groups::RowId;
use fs_assimilate::groups::jacobian::{JacobianBatch, JacobianRecord, SharedNoiseModel};
use fs_assimilate::groups::stream::SharedNoiseStream;
use fs_assimilate::{AssimError, Belief, Observation, assimilate};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};

fn with_cx<R>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 6401,
                kernel_id: 64,
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

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-12 * (1.0 + a.abs().max(b.abs())),
        "{a} != {b}"
    );
}

fn zero_noise(n: usize, cx: &Cx<'_>) -> FactorBelief {
    FactorBelief::diagonal(vec![0.0; n], &vec![0.0; n], cx).unwrap()
}

#[test]
fn square_root_forecast_matches_full_covariance_formula() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        for n in 1..=8 {
            let p: Vec<Vec<_>> = (0..n)
                .map(|i| {
                    (0..n)
                        .map(|j| if i == j { 2.0 + i as f64 } else { 0.125 })
                        .collect()
                })
                .collect();
            let q: Vec<Vec<_>> = (0..n)
                .map(|i| (0..n).map(|j| if i == j { 0.25 } else { -0.01 }).collect())
                .collect();
            let m: Vec<_> = (0..n).map(|i| i as f64 * 0.5 - 1.0).collect();
            let f: Vec<Vec<_>> = (0..n)
                .map(|i| {
                    (0..n)
                        .map(|j| {
                            if i == j {
                                0.8
                            } else {
                                ((i * 7 + j * 3) % 5) as f64 * 0.1 - 0.2
                            }
                        })
                        .collect()
                })
                .collect();
            let b: Vec<_> = (0..n).map(|i| -0.3 + i as f64 * 0.2).collect();
            let prior =
                FactorBelief::from_belief(&Belief::new(m.clone(), p.clone(), cx).unwrap(), cx)
                    .unwrap();
            let noise =
                FactorBelief::from_belief(&Belief::new(vec![0.0; n], q.clone(), cx).unwrap(), cx)
                    .unwrap();
            let dynamics =
                LinearDynamics::new(f.clone(), b.clone(), noise, "linear-test", cx).unwrap();
            let result = dynamics.predict(&prior, cx).unwrap();
            assert_eq!(result, dynamics.predict(&prior, cx).unwrap());
            let covariance = result.to_dense_covariance();
            for i in 0..n {
                close(
                    result.mean()[i],
                    b[i] + (0..n).map(|k| f[i][k] * m[k]).sum::<f64>(),
                );
                for j in 0..n {
                    let mut expected = q[i][j];
                    for k in 0..n {
                        for l in 0..n {
                            expected += f[i][k] * p[k][l] * f[j][l];
                        }
                    }
                    close(covariance[i][j], expected);
                }
            }
        }
    });
}

#[test]
fn repeated_prediction_and_measurement_follow_scalar_kalman_recursion() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let dynamics = LinearDynamics::new(
            vec![vec![1.2]],
            vec![0.3],
            FactorBelief::diagonal(vec![0.0], &[0.4], cx).unwrap(),
            "scalar-model",
            cx,
        )
        .unwrap();
        let mut factor = FactorBelief::diagonal(vec![-0.2], &[2.0], cx).unwrap();
        let (mut mean, mut variance) = (-0.2, 2.0);
        for i in 0..20 {
            let y = i as f64 * 0.1;
            mean = 1.2 * mean + 0.3;
            variance = 1.44 * variance + 0.4;
            let gain = variance / (variance + 0.2);
            mean += gain * (y - mean);
            variance *= 0.2 / (variance + 0.2);
            let predicted = dynamics.predict(&factor, cx).unwrap();
            factor = assimilate_scalar(
                &predicted,
                &Observation::new(vec![1.0], y, 0.2, "sensor").unwrap(),
                cx,
            )
            .unwrap()
            .belief()
            .clone();
            close(factor.mean()[0], mean);
            close(factor.variance(0).unwrap(), variance);
        }
    });
}

#[test]
fn exact_null_rows_do_not_discard_remaining_variance() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let prior = FactorBelief::diagonal(vec![1.0, 2.0, 3.0], &[2.0, 4.0, 8.0], cx).unwrap();
        for f in [
            vec![
                vec![1.0, 0.0, 0.0],
                vec![0.0, 0.0, 0.0],
                vec![0.0, 0.0, 0.0],
            ],
            vec![
                vec![0.0, 0.0, 1.0],
                vec![0.0, 0.0, 0.0],
                vec![0.0, 1.0, 0.0],
            ],
            vec![
                vec![0.0, 0.0, 0.0],
                vec![0.0, 0.0, 0.0],
                vec![0.0, 0.0, 0.0],
            ],
        ] {
            let dynamics = LinearDynamics::new(
                f.clone(),
                vec![0.0; 3],
                zero_noise(3, cx),
                "rank-deficient",
                cx,
            )
            .unwrap();
            let predicted = dynamics.predict(&prior, cx).unwrap();
            for i in 0..3 {
                let expected = (0..3)
                    .map(|j| f[i][j] * f[i][j] * prior.diag(j).unwrap())
                    .sum::<f64>();
                close(predicted.variance(i).unwrap(), expected);
                if expected == 0.0 {
                    assert_eq!(predicted.variance(i), Some(0.0));
                }
            }
        }
    });
}

#[test]
fn rescaled_roots_avoid_squaring_transition_before_scaling_covariance() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let prior = FactorBelief::diagonal(vec![0.0; 2], &[1e-300, 1e300], cx).unwrap();
        let dynamics = LinearDynamics::new(
            vec![vec![1e150, 0.0], vec![0.0, 1e-150]],
            vec![0.0; 2],
            zero_noise(2, cx),
            "rescaled",
            cx,
        )
        .unwrap();
        let prediction = dynamics.predict(&prior, cx).unwrap();
        close(prediction.variance(0).unwrap(), 1.0);
        close(prediction.variance(1).unwrap(), 1.0);
        assert_eq!(prediction.to_dense_covariance()[0][1], 0.0);
    });
}

fn sample(model: &SharedNoiseModel, id: &str, value: f64, cx: &Cx<'_>) -> JacobianBatch {
    let record = JacobianRecord::new(
        Observation::new(vec![1.0], value, 0.25, "sensor").unwrap(),
        RowId::try_new(id).unwrap(),
        model,
        vec![1.0],
    )
    .unwrap();
    JacobianBatch::new(model.clone(), vec![record], cx).unwrap()
}

#[test]
fn dynamic_stream_retains_calibration_cross_covariance_and_replay_protection() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let model =
            SharedNoiseModel::new("calibration", vec!["offset".into()], vec![vec![1.0]], cx)
                .unwrap();
        let mut stream =
            SharedNoiseStream::new(&Belief::scalar(1.0, 4.0).unwrap(), &model, 4, cx).unwrap();
        let first = sample(&model, "first", 3.0, cx);
        stream.assimilate(&first, cx).unwrap();
        let before = stream.clone();
        let p = before.joint().to_dense_covariance();
        assert!(p[0][1].abs() > 0.1);
        let dynamics = LinearDynamics::new(
            vec![vec![2.0]],
            vec![0.5],
            FactorBelief::diagonal(vec![0.0], &[0.2], cx).unwrap(),
            "growth",
            cx,
        )
        .unwrap();
        stream.forecast_physical(&dynamics, 0.0, 0.25, cx).unwrap();
        assert_eq!(stream.time(), Some(0.25));
        assert_eq!(stream.accepted_rows(), 1);
        assert_ne!(stream.history_id(), before.history_id());
        let expected = Belief::new(
            vec![
                2.0 * before.joint().mean()[0] + 0.5,
                before.joint().mean()[1],
            ],
            vec![
                vec![4.0 * p[0][0] + 0.2, 2.0 * p[0][1]],
                vec![2.0 * p[1][0], p[1][1]],
            ],
            cx,
        )
        .unwrap();
        let actual = stream.joint().to_dense_covariance();
        assert_eq!(stream.joint().mean()[1], before.joint().mean()[1]);
        assert_eq!(stream.joint().diag(1), before.joint().diag(1));
        for i in 0..2 {
            for j in 0..2 {
                close(actual[i][j], expected.covariance()[i][j]);
            }
        }
        let predicted = stream.clone();
        stream.assimilate(&first, cx).unwrap();
        assert_eq!(stream, predicted);
        stream
            .assimilate(&sample(&model, "second", 4.0, cx), cx)
            .unwrap();
        let oracle = assimilate(
            &expected,
            &Observation::new(vec![1.0, 1.0], 4.0, 0.25, "joint").unwrap(),
            cx,
        )
        .unwrap();
        let actual = stream.joint().to_dense_covariance();
        for i in 0..2 {
            close(stream.joint().mean()[i], oracle.mean()[i]);
            for j in 0..2 {
                close(actual[i][j], oracle.covariance()[i][j]);
            }
        }
        let mut replay = before;
        replay.forecast_physical(&dynamics, 0.0, 0.25, cx).unwrap();
        assert_eq!(replay, predicted);
    });
}

#[test]
fn invalid_dynamics_and_time_frontiers_refuse_without_mutation() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        assert!(
            LinearDynamics::new(
                vec![vec![1.0, 0.0]],
                vec![0.0],
                zero_noise(1, cx),
                "bad",
                cx
            )
            .is_err()
        );
        assert!(
            LinearDynamics::new(
                vec![vec![f64::NAN]],
                vec![0.0],
                zero_noise(1, cx),
                "bad",
                cx
            )
            .is_err()
        );
        assert!(
            LinearDynamics::new(
                vec![vec![1.0]],
                vec![f64::INFINITY],
                zero_noise(1, cx),
                "bad",
                cx
            )
            .is_err()
        );
        let shifted_noise = FactorBelief::diagonal(vec![1.0], &[0.2], cx).unwrap();
        assert!(LinearDynamics::new(vec![vec![1.0]], vec![0.0], shifted_noise, "bad", cx).is_err());
        let overflowing_covariance = LinearDynamics::new(
            vec![vec![1e200, 0.0], vec![1.0, 0.0]],
            vec![0.0; 2],
            zero_noise(2, cx),
            "finite-pivots-overflow",
            cx,
        )
        .unwrap();
        let prior = FactorBelief::diagonal(vec![0.0; 2], &[1.0; 2], cx).unwrap();
        assert!(matches!(
            overflowing_covariance.predict(&prior, cx),
            Err(AssimError::NonFiniteComputation {
                stage: "forecast marginal variance"
            })
        ));
        let dynamics = LinearDynamics::new(
            vec![vec![1.0]],
            vec![0.0],
            zero_noise(1, cx),
            "identity",
            cx,
        )
        .unwrap();
        let model =
            SharedNoiseModel::new("calibration", vec!["offset".into()], vec![vec![1.0]], cx)
                .unwrap();
        let mut stream =
            SharedNoiseStream::new(&Belief::scalar(0.0, 1.0).unwrap(), &model, 2, cx).unwrap();
        stream.forecast_physical(&dynamics, 0.0, 1.0, cx).unwrap();
        let before = stream.clone();
        for (from, to) in [
            (0.0, 1.0),
            (0.5, 2.0),
            (1.0, 1.0),
            (1.0, 0.0),
            (f64::NAN, 2.0),
            (1.0, f64::INFINITY),
        ] {
            assert!(matches!(
                stream.forecast_physical(&dynamics, from, to, cx),
                Err(AssimError::InvalidForecast { .. })
            ));
            assert_eq!(stream, before);
        }
        stream.forecast_physical(&dynamics, 1.0, 2.0, cx).unwrap();
        assert_eq!(stream.time(), Some(2.0));
    });
}

#[test]
fn forecast_cancellation_quota_and_overflow_are_atomic_and_retryable() {
    let (source, dynamics, prior) = with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let prior = FactorBelief::diagonal(vec![0.5], &[1.0], cx).unwrap();
        let dynamics =
            LinearDynamics::new(vec![vec![1.1]], vec![0.1], zero_noise(1, cx), "step", cx).unwrap();
        let model =
            SharedNoiseModel::new("calibration", vec!["offset".into()], vec![vec![1.0]], cx)
                .unwrap();
        (
            SharedNoiseStream::new(&Belief::scalar(0.5, 1.0).unwrap(), &model, 2, cx).unwrap(),
            dynamics,
            prior,
        )
    });
    for budget in [
        Budget::INFINITE.with_poll_quota(0),
        Budget::INFINITE.with_poll_quota(1),
        Budget::INFINITE.with_cost_quota(1),
    ] {
        with_cx(&CancelGate::new(), budget, |cx| {
            assert!(dynamics.predict(&prior, cx).is_err());
            let mut candidate = source.clone();
            assert!(
                candidate
                    .forecast_physical(&dynamics, 0.0, 1.0, cx)
                    .is_err()
            );
            assert_eq!(candidate, source);
        });
    }
    let gate = CancelGate::new();
    gate.request();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let mut candidate = source.clone();
        assert!(
            candidate
                .forecast_physical(&dynamics, 0.0, 1.0, cx)
                .is_err()
        );
        assert_eq!(candidate, source);
    });
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let mut candidate = source.clone();
        let overflow = LinearDynamics::new(
            vec![vec![1e200]],
            vec![0.0],
            zero_noise(1, cx),
            "overflow",
            cx,
        )
        .unwrap();
        assert!(
            candidate
                .forecast_physical(&overflow, 0.0, 1.0, cx)
                .is_err()
        );
        assert_eq!(candidate, source);
        candidate
            .forecast_physical(&dynamics, 0.0, 1.0, cx)
            .unwrap();
        assert_eq!(candidate.time(), Some(1.0));
    });
}

#[test]
fn mixed_physical_dynamics_keep_the_full_learned_nuisance_block() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let model = SharedNoiseModel::new(
            "two-calibrations",
            vec!["offset".into(), "tilt".into()],
            vec![vec![0.7, 0.1], vec![0.1, 0.4]],
            cx,
        )
        .unwrap();
        let prior = Belief::new(vec![0.4, -0.2], vec![vec![2.0, 0.3], vec![0.3, 1.0]], cx).unwrap();
        let mut stream = SharedNoiseStream::new(&prior, &model, 4, cx).unwrap();
        let record = JacobianRecord::new(
            Observation::new(vec![1.0, 0.5], 2.0, 0.2, "sensor").unwrap(),
            RowId::try_new("first").unwrap(),
            &model,
            vec![0.3, -1.2],
        )
        .unwrap();
        stream
            .assimilate(&JacobianBatch::new(model, vec![record], cx).unwrap(), cx)
            .unwrap();
        let before = stream.clone();
        let old_p = before.joint().to_dense_covariance();
        let f = vec![vec![1.0, 0.2], vec![-0.1, 0.9]];
        let q = vec![vec![0.2, 0.05], vec![0.05, 0.3]];
        let noise =
            FactorBelief::from_belief(&Belief::new(vec![0.0; 2], q.clone(), cx).unwrap(), cx)
                .unwrap();
        let dynamics =
            LinearDynamics::new(f.clone(), vec![0.1, -0.2], noise, "coupled-motion", cx).unwrap();
        stream.forecast_physical(&dynamics, 10.0, 11.0, cx).unwrap();
        let a = [
            [1.0, 0.2, 0.0, 0.0],
            [-0.1, 0.9, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let actual = stream.joint().to_dense_covariance();
        for i in 0..4 {
            for j in 0..4 {
                let mut expected = if i < 2 && j < 2 { q[i][j] } else { 0.0 };
                for k in 0..4 {
                    for l in 0..4 {
                        expected += a[i][k] * old_p[k][l] * a[j][l];
                    }
                }
                close(actual[i][j], expected);
                if i >= 2 && j >= 2 {
                    assert_eq!(actual[i][j], old_p[i][j]);
                }
            }
        }
    });
}

#[test]
fn forecast_history_binds_models_intervals_and_invalidates_sensor_plans() {
    use fs_assimilate::groups::stream::design::SharedNoiseProbe;
    use fs_assimilate::groups::stream::design::batch::{CostedProbe, ProbeBatchBudget};
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let model =
            SharedNoiseModel::new("calibration", vec!["offset".into()], vec![vec![1.0]], cx)
                .unwrap();
        let original =
            SharedNoiseStream::new(&Belief::scalar(0.0, 1.0).unwrap(), &model, 4, cx).unwrap();
        let probe = CostedProbe::new(
            SharedNoiseProbe::new(
                &model,
                RowId::try_new("next").unwrap(),
                vec![1.0],
                vec![1.0],
                0.25,
            )
            .unwrap(),
            1,
        )
        .unwrap();
        let plan = original
            .plan_probe_batch(&[1.0], &[probe], ProbeBatchBudget::new(1, 1).unwrap(), cx)
            .unwrap();
        plan.check_current(&original).unwrap();
        let a = LinearDynamics::new(vec![vec![1.0]], vec![0.0], zero_noise(1, cx), "model-a", cx)
            .unwrap();
        let b = LinearDynamics::new(vec![vec![1.0]], vec![0.0], zero_noise(1, cx), "model-b", cx)
            .unwrap();
        assert_ne!(a.content_id(), b.content_id());
        let mut first = original.clone();
        first.forecast_physical(&a, 0.0, 1.0, cx).unwrap();
        assert!(plan.check_current(&first).is_err());
        let mut second = original.clone();
        second.forecast_physical(&b, 0.0, 1.0, cx).unwrap();
        assert_eq!(first.joint(), second.joint());
        assert_ne!(first.history_id(), second.history_id());
        let mut later = original;
        later.forecast_physical(&a, 0.0, 2.0, cx).unwrap();
        assert_eq!(first.joint(), later.joint());
        assert_ne!(first.history_id(), later.history_id());
    });
}


#[test]
fn nonlinear_vector_prediction_matches_analytic_jacobian_and_affine_limit() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let covariance = vec![vec![2.0, 0.5], vec![0.5, 1.0]];
        let noise_covariance = vec![vec![0.125, 0.015625], vec![0.015625, 0.25]];
        let prior = FactorBelief::from_belief(
            &Belief::new(vec![2.0, -1.0], covariance.clone(), cx).unwrap(), cx,
        ).unwrap();
        let noise = FactorBelief::from_belief(
            &Belief::new(vec![0.0; 2], noise_covariance.clone(), cx).unwrap(), cx,
        ).unwrap();
        let dynamics = NonlinearDynamics::new(
            noise.clone(), vec![0.25, 0.125], "quadratic-transition", cx,
        ).unwrap();
        let mut calls = 0;
        let predicted = dynamics.predict_fd(&prior, |x, _| {
            calls += 1;
            Ok(vec![x[0] * x[0] + 0.5 * x[1], x[0] * x[1]])
        }, cx).unwrap();
        assert_eq!(calls, 5, "one baseline and two vector probes per component");
        assert_eq!(predicted.mean(), &[3.5, -2.0]);
        let jacobian = [[4.0, 0.5], [-1.0, 2.0]];
        let actual = predicted.to_dense_covariance();
        for i in 0..2 {
            for j in 0..2 {
                let mut expected = noise_covariance[i][j];
                for k in 0..2 {
                    for l in 0..2 {
                        expected += jacobian[i][k] * covariance[k][l] * jacobian[j][l];
                    }
                }
                close(actual[i][j], expected);
            }
        }
        let affine = NonlinearDynamics::new(
            noise.clone(), vec![0.25, 0.125], "affine-oracle", cx,
        ).unwrap();
        let linear = LinearDynamics::new(
            vec![vec![0.5, 0.25], vec![-1.0, 2.0]], vec![0.125, -0.5],
            noise, "affine-oracle", cx,
        ).unwrap();
        let local = affine.predict_fd(&prior, |x, _| {
            Ok(vec![0.5 * x[0] + 0.25 * x[1] + 0.125, -x[0] + 2.0 * x[1] - 0.5])
        }, cx).unwrap();
        assert_eq!(local, linear.predict(&prior, cx).unwrap());
    });
}

#[test]
fn nonlinear_mean_is_not_lost_by_constructing_a_large_affine_offset() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let prior = FactorBelief::diagonal(vec![1.0], &[1e-32], cx).unwrap();
        let dynamics = NonlinearDynamics::new(
            zero_noise(1, cx), vec![1.0 / 1024.0], "locally-steep", cx,
        ).unwrap();
        let predicted = dynamics.predict_fd(
            &prior, |x, _| Ok(vec![1.0 + 1e16 * (x[0] - 1.0)]), cx,
        ).unwrap();
        assert_eq!(predicted.mean(), &[1.0]);
        close(predicted.variance(0).unwrap(), 1.0);
    });
}

#[test]
fn nonlinear_stream_keeps_learned_calibration_and_matches_next_joint_update() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let model = SharedNoiseModel::new(
            "nonlinear-calibration", vec!["offset".into()], vec![vec![1.0]], cx,
        ).unwrap();
        let mut stream = SharedNoiseStream::new(
            &Belief::scalar(1.0, 4.0).unwrap(), &model, 4, cx,
        ).unwrap();
        let first = sample(&model, "nonlinear-first", 3.0, cx);
        stream.assimilate(&first, cx).unwrap();
        let before = stream.clone();
        let p = before.joint().to_dense_covariance();
        let m = before.joint().mean()[0];
        let derivative = 2.0 * m;
        let dynamics = NonlinearDynamics::new(
            FactorBelief::diagonal(vec![0.0], &[0.2], cx).unwrap(),
            vec![0.25], "square-dynamics", cx,
        ).unwrap();
        stream.forecast_physical_nonlinear(
            &dynamics, 0.0, 1.0, |x, _| Ok(vec![x[0] * x[0]]), cx,
        ).unwrap();
        let expected = Belief::new(
            vec![m * m, before.joint().mean()[1]],
            vec![
                vec![derivative * derivative * p[0][0] + 0.2, derivative * p[0][1]],
                vec![derivative * p[1][0], p[1][1]],
            ], cx,
        ).unwrap();
        let predicted = stream.clone();
        assert_eq!(stream.time(), Some(1.0));
        assert_ne!(stream.history_id(), before.history_id());
        assert_eq!(stream.joint().diag(1), before.joint().diag(1));
        for i in 0..2 {
            close(stream.joint().mean()[i], expected.mean()[i]);
            for j in 0..2 {
                close(stream.joint().to_dense_covariance()[i][j], expected.covariance()[i][j]);
            }
        }
        stream.assimilate(&first, cx).unwrap();
        assert_eq!(stream, predicted, "time evolution must not reissue an old reading");
        let mut calls = 0;
        assert!(stream.forecast_physical_nonlinear(&dynamics, 0.0, 1.0, |x, _| {
            calls += 1;
            Ok(x.to_vec())
        }, cx).is_err());
        assert_eq!(calls, 0, "a bad frontier is refused before model work");
        assert_eq!(stream, predicted);
        let next = sample(&model, "nonlinear-second", 6.0, cx);
        stream.assimilate(&next, cx).unwrap();
        let oracle = assimilate(
            &expected, &Observation::new(vec![1.0, 1.0], 6.0, 0.25, "joint").unwrap(), cx,
        ).unwrap();
        let covariance = stream.joint().to_dense_covariance();
        for i in 0..2 {
            close(stream.joint().mean()[i], oracle.mean()[i]);
            for j in 0..2 { close(covariance[i][j], oracle.covariance()[i][j]); }
        }
    });
}

#[test]
fn nonlinear_probe_failures_and_cancellation_never_advance_the_stream() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let model = SharedNoiseModel::new(
            "probe-calibration", vec!["offset".into()], vec![vec![1.0]], cx,
        ).unwrap();
        let mut stream = SharedNoiseStream::new(
            &Belief::scalar(1.0, 1.0).unwrap(), &model, 4, cx,
        ).unwrap();
        let dynamics = NonlinearDynamics::new(
            zero_noise(1, cx), vec![0.25], "domain-limited", cx,
        ).unwrap();
        let saved = stream.clone();
        for failure in 0..3 {
            let result = stream.forecast_physical_nonlinear(&dynamics, 0.0, 1.0, |x, _| {
                if x[0] < 1.0 {
                    match failure {
                        0 => Err(AssimError::InvalidForecast { reason: "outside model domain" }),
                        1 => Ok(vec![f64::NAN]),
                        _ => Ok(vec![]),
                    }
                } else { Ok(vec![x[0] * x[0]]) }
            }, cx);
            assert!(result.is_err());
            assert_eq!(stream, saved);
        }
        let gate = CancelGate::new();
        with_cx(&gate, Budget::INFINITE, |cancelled| {
            assert!(matches!(stream.forecast_physical_nonlinear(
                &dynamics, 0.0, 1.0, |x, _| { gate.request(); Ok(x.to_vec()) }, cancelled,
            ), Err(AssimError::Cancelled { .. })));
        });
        assert_eq!(stream, saved);
        let mut calls = 0;
        with_cx(&CancelGate::new(), Budget { cost_quota: Some(1), ..Budget::INFINITE }, |tiny| {
            assert!(stream.forecast_physical_nonlinear(&dynamics, 0.0, 1.0, |x, _| {
                calls += 1;
                Ok(x.to_vec())
            }, tiny).is_err());
        });
        assert_eq!(calls, 0);
        assert_eq!(stream, saved);
        let large = FactorBelief::diagonal(vec![1e30], &[1.0], cx).unwrap();
        assert!(matches!(dynamics.predict_fd(&large, |x, _| {
            calls += 1;
            Ok(x.to_vec())
        }, cx), Err(AssimError::FiniteDifferenceStepUnrepresentable { component: 0 })));
        assert_eq!(calls, 0);
        stream.forecast_physical_nonlinear(
            &dynamics, 0.0, 1.0, |x, _| Ok(vec![x[0] * x[0]]), cx,
        ).unwrap();
        assert_eq!(stream.time(), Some(1.0), "a failed attempt leaves the interval retryable");
    });
}

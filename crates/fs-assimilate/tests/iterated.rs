//! G1/G3/G4: nonlinear iteration against linear fusion and a scalar MAP oracle.

use std::cell::Cell;

use fs_assimilate::nonlinear::iterated::{
    IndependentMeasurements, IteratedOutcome, IteratedSettings, ModelEvaluation,
    assimilate_iterated,
};
use fs_assimilate::{AssimError, Belief, Observation, assimilate_all};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};

fn with_cx<R>(gate: &CancelGate, run: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 42,
                kernel_id: 12,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        )
        .with_time_source(&clock);
        run(&cx)
    });
    assert!(pool.stats().quiescent());
    result
}

fn settings() -> IteratedSettings {
    IteratedSettings {
        max_iterations: 64,
        max_model_evaluations: 256,
        max_backtracks: 16,
        gradient_tolerance: 1e-7,
        step_tolerance: 1e-14,
    }
}

fn measurement(value: f64, variance: f64) -> IndependentMeasurements {
    IndependentMeasurements::new(
        vec![value],
        vec![variance],
        vec!["sensor".into()],
        "square-v1",
    )
    .expect("valid independent measurement")
}

fn square(state: &[f64], _: &Cx<'_>) -> Result<ModelEvaluation, AssimError> {
    Ok(ModelEvaluation::Value {
        predictions: vec![state[0] * state[0]],
        jacobian: vec![vec![2.0 * state[0]]],
    })
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

fn belief_close(actual: &Belief, expected: &Belief) {
    for (&a, &b) in actual.mean().iter().zip(expected.mean()) {
        close(a, b, 1e-11);
    }
    for (a, b) in actual.covariance().iter().zip(expected.covariance()) {
        for (&a, &b) in a.iter().zip(b) {
            close(a, b, 1e-11);
        }
    }
}

#[test]
fn affine_multirow_matches_one_linear_batch_without_reusing_evidence() {
    with_cx(&CancelGate::new(), |cx| {
        let prior = Belief::new(vec![0.5, -1.0], vec![vec![2.0, 0.5], vec![0.5, 1.0]], cx)
            .expect("correlated positive definite prior");
        let inputs = IndependentMeasurements::new(
            vec![2.0, -1.0],
            vec![0.5, 2.0],
            vec!["sensor-a".into(), "sensor-b".into()],
            "affine-v1",
        )
        .unwrap();
        let linear = [
            Observation::new(vec![1.0, 2.0], 2.0 - 0.25, 0.5, "sensor-a").unwrap(),
            Observation::new(vec![-2.0, 1.0], -1.0 + 0.5, 2.0, "sensor-b").unwrap(),
        ];
        let oracle = assimilate_all(&prior, &linear, cx).unwrap();
        for max_iterations in [2, 64] {
            let mut policy = settings();
            policy.max_iterations = max_iterations;
            let model = |x: &[f64], _: &Cx<'_>| {
                Ok(ModelEvaluation::Value {
                    predictions: vec![x[0] + 2.0 * x[1] + 0.25, -2.0 * x[0] + x[1] - 0.5],
                    jacobian: vec![vec![1.0, 2.0], vec![-2.0, 1.0]],
                })
            };
            let result = assimilate_iterated(&prior, &inputs, policy, model, cx).unwrap();
            assert_eq!(result.outcome(), IteratedOutcome::Stationary);
            assert_eq!(result.trace().len(), 2);
            assert_eq!(result.model_evaluations(), 2);
            belief_close(result.posterior().unwrap(), &oracle);
        }
    });
}

#[test]
fn square_observation_reaches_independently_solved_prior_regularized_optimum() {
    // J(x) = (x-1)^2/2 + (4-x^2)^2/(2*0.25). Its derivative is
    // strictly increasing on [1.5,2], and changes sign across that bracket.
    let (mut lo, mut hi) = (1.5, 2.0);
    for _ in 0..80 {
        let x = 0.5 * (lo + hi);
        if (x - 1.0) + 8.0 * x * (x * x - 4.0) > 0.0 {
            hi = x;
        } else {
            lo = x;
        }
    }
    let optimum = 0.5 * (lo + hi);
    with_cx(&CancelGate::new(), |cx| {
        let prior = Belief::scalar(1.0, 1.0).unwrap();
        let result =
            assimilate_iterated(&prior, &measurement(4.0, 0.25), settings(), square, cx).unwrap();
        assert_eq!(result.outcome(), IteratedOutcome::Stationary);
        close(result.state()[0], optimum, 1e-8);
        let trace = result.trace();
        assert!(trace.len() > 2, "must relinearize the nonlinear model");
        assert!(trace.windows(2).all(|w| w[1].objective < w[0].objective));
        assert!(result.trace().last().unwrap().gradient_norm <= settings().gradient_tolerance);
        // This is the local Gauss-Newton covariance, not the exact posterior
        // variance or the full Hessian inverse of the non-Gaussian objective.
        let local_variance = 1.0 / (1.0 + 16.0 * optimum * optimum);
        close(
            result.posterior().unwrap().variance(0).unwrap(),
            local_variance,
            1e-9,
        );
    });
}

#[test]
fn domain_backtracking_shortens_trials_and_retains_only_accepted_states() {
    with_cx(&CancelGate::new(), |cx| {
        let prior = Belief::scalar(0.1, 1.0).unwrap();
        let calls = Cell::new(0);
        let outside = Cell::new(0);
        let model = |x: &[f64], cx: &Cx<'_>| {
            calls.set(calls.get() + 1);
            if !(0.0..=3.0).contains(&x[0]) {
                outside.set(outside.get() + 1);
                Ok(ModelEvaluation::OutsideDomain)
            } else {
                square(x, cx)
            }
        };
        let inputs = measurement(4.0, 0.01);
        let result = assimilate_iterated(&prior, &inputs, settings(), model, cx).unwrap();
        assert_eq!(result.outcome(), IteratedOutcome::Stationary);
        assert!(outside.get() >= 3);
        assert!(result.rejected_trials() >= outside.get());
        assert_eq!(result.model_evaluations(), calls.get());
        assert!(result.trace()[1].step_fraction < 1.0);
        let trace = result.trace();
        assert!(trace.windows(2).all(|w| w[1].objective < w[0].objective));
        assert!((1.99..2.0).contains(&result.state()[0]));
    });
}

#[test]
fn singular_correlated_prior_preserves_its_support() {
    with_cx(&CancelGate::new(), |cx| {
        let prior = Belief::new(vec![0.0, 0.0], vec![vec![1.0, 1.0], vec![1.0, 1.0]], cx).unwrap();
        let inputs = measurement(2.0, 1.0);
        let model = |x: &[f64], _: &Cx<'_>| {
            Ok(ModelEvaluation::Value {
                predictions: vec![x[0] + x[1]],
                jacobian: vec![vec![1.0, 1.0]],
            })
        };
        let result = assimilate_iterated(&prior, &inputs, settings(), model, cx).unwrap();
        assert_eq!(result.outcome(), IteratedOutcome::Stationary);
        let posterior = result.posterior().unwrap();
        for &x in posterior.mean() {
            close(x, 0.8, 1e-12);
        }
        for row in posterior.covariance() {
            for &value in row {
                close(value, 0.2, 1e-12);
            }
        }
        close(posterior.mean()[0] - posterior.mean()[1], 0.0, 1e-14);
    });
}

#[test]
fn numerical_caps_and_exhausted_domain_search_publish_no_posterior() {
    with_cx(&CancelGate::new(), |cx| {
        let prior = Belief::scalar(1.0, 1.0).unwrap();
        let inputs = measurement(4.0, 0.25);
        for (iterations, calls, expected, consumed) in [
            (0, 256, IteratedOutcome::IterationLimit, 1),
            (64, 1, IteratedOutcome::EvaluationLimit, 1),
            (1, 256, IteratedOutcome::IterationLimit, 2),
        ] {
            let mut policy = settings();
            policy.max_iterations = iterations;
            policy.max_model_evaluations = calls;
            let result = assimilate_iterated(&prior, &inputs, policy, square, cx).unwrap();
            assert_eq!(result.outcome(), expected);
            assert_eq!(result.model_evaluations(), consumed);
            assert!(result.posterior().is_none());
        }
        for (calls, expected, rejected) in [
            (256, IteratedOutcome::LineSearchStalled, 3),
            (2, IteratedOutcome::EvaluationLimit, 1),
        ] {
            let mut policy = settings();
            policy.max_model_evaluations = calls;
            policy.max_backtracks = 2;
            let result = assimilate_iterated(
                &prior,
                &inputs,
                policy,
                |x, cx| {
                    if x == prior.mean() {
                        square(x, cx)
                    } else {
                        Ok(ModelEvaluation::OutsideDomain)
                    }
                },
                cx,
            )
            .unwrap();
            assert_eq!(result.outcome(), expected);
            assert!(result.posterior().is_none());
            assert_eq!(result.state(), prior.mean());
            assert_eq!(result.trace().len(), 1);
            assert_eq!(result.rejected_trials(), rejected);
        }
    });
}

#[test]
fn cancellation_before_and_after_callbacks_is_atomic() {
    for pre_cancelled in [true, false] {
        let gate = CancelGate::new();
        let calls = Cell::new(0);
        if pre_cancelled {
            gate.request();
        }
        with_cx(&gate, |cx| {
            let prior = Belief::scalar(1.0, 1.0).unwrap();
            let original = prior.clone();
            let model = |x: &[f64], cx: &Cx<'_>| {
                calls.set(calls.get() + 1);
                gate.request();
                square(x, cx)
            };
            let inputs = measurement(4.0, 0.25);
            let result = assimilate_iterated(&prior, &inputs, settings(), model, cx);
            assert!(matches!(result, Err(AssimError::Cancelled { .. })));
            assert_eq!(prior, original);
            assert_eq!(calls.get(), usize::from(!pre_cancelled));
        });
    }
}

#[test]
fn malformed_measurements_models_and_policies_are_refused() {
    for (values, variances, instruments) in [
        (vec![], vec![], vec![]),
        (vec![1.0], vec![], vec!["sensor".into()]),
        (vec![1.0], vec![0.0], vec!["sensor".into()]),
        (vec![f64::NAN], vec![1.0], vec!["sensor".into()]),
        (
            vec![1.0, 2.0],
            vec![1.0, 1.0],
            vec!["sensor".into(), "sensor".into()],
        ),
    ] {
        assert!(IndependentMeasurements::new(values, variances, instruments, "model-v1").is_err());
    }
    with_cx(&CancelGate::new(), |cx| {
        let prior = Belief::scalar(1.0, 1.0).unwrap();
        let inputs = measurement(4.0, 0.25);
        for (predictions, jacobian) in [
            (vec![], vec![vec![1.0]]),
            (vec![1.0], vec![]),
            (vec![1.0], vec![vec![1.0, 2.0]]),
            (vec![1.0], vec![vec![f64::NAN]]),
            (vec![f64::INFINITY], vec![vec![1.0]]),
        ] {
            let result = assimilate_iterated(
                &prior,
                &inputs,
                settings(),
                |_, _| {
                    Ok(ModelEvaluation::Value {
                        predictions: predictions.clone(),
                        jacobian: jacobian.clone(),
                    })
                },
                cx,
            );
            assert!(result.is_err());
        }
        let outside = |_: &[f64], _: &Cx<'_>| Ok(ModelEvaluation::OutsideDomain);
        assert!(assimilate_iterated(&prior, &inputs, settings(), outside, cx).is_err());
        let calls = Cell::new(0);
        let mut policy = settings();
        policy.max_model_evaluations = 0;
        let result = assimilate_iterated(
            &prior,
            &inputs,
            policy,
            |x, cx| {
                calls.set(calls.get() + 1);
                square(x, cx)
            },
            cx,
        );
        assert!(result.is_err());
        assert_eq!(calls.get(), 0);
    });
}

#[test]
fn deterministic_replay_retains_state_covariance_and_receipt() {
    let run = || {
        with_cx(&CancelGate::new(), |cx| {
            assimilate_iterated(
                &Belief::scalar(1.0, 1.0).unwrap(),
                &measurement(4.0, 0.25),
                settings(),
                square,
                cx,
            )
            .unwrap()
        })
    };
    let first = run();
    assert_eq!(first.outcome(), IteratedOutcome::Stationary);
    assert_eq!(first, run());
    assert!(!first.receipt_id().is_empty());
}

//! G1/G3/G4: correlated nonlinear likelihoods, shared-noise floors and refusal.

use std::cell::Cell;

use fs_assimilate::nonlinear::iterated::{
    CorrelatedMeasurements, IndependentMeasurements, IteratedOutcome, IteratedSettings,
    ModelEvaluation, assimilate_correlated_iterated, assimilate_iterated,
};
use fs_assimilate::{
    AssimError, Belief, Observation, ObservationBatch, RobustObservation,
    assimilate_observation_batch,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};

fn with_cx<R>(gate: &CancelGate, run: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| {
        let key = StreamKey {
            seed: 73,
            kernel_id: 19,
            tile: 0,
            iteration: 0,
        };
        let cx = Cx::new(gate, arena, key, Budget::INFINITE, ExecMode::Deterministic)
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

fn instruments(n: usize) -> Vec<String> {
    (0..n).map(|i| format!("sensor-{i}")).collect()
}

fn inputs(values: Vec<f64>, covariance: Vec<Vec<f64>>, cx: &Cx<'_>) -> CorrelatedMeasurements {
    let ids = instruments(values.len());
    CorrelatedMeasurements::new(values, covariance, ids, "calibration-v1", cx).unwrap()
}

fn noise() -> Vec<Vec<f64>> {
    vec![vec![0.25, 0.125], vec![0.125, 0.5]]
}

fn nonlinear(x: &[f64], _: &Cx<'_>) -> Result<ModelEvaluation, AssimError> {
    Ok(ModelEvaluation::Value {
        predictions: vec![x[0] * x[0], x[0]],
        jacobian: vec![vec![2.0 * x[0]], vec![1.0]],
    })
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

fn belief_close(actual: &Belief, expected: &Belief, tolerance: f64) {
    assert_eq!(actual.dim(), expected.dim());
    for (&a, &b) in actual.mean().iter().zip(expected.mean()) {
        close(a, b, tolerance);
    }
    for (a, b) in actual.covariance().iter().zip(expected.covariance()) {
        for (&a, &b) in a.iter().zip(b) {
            close(a, b, tolerance);
        }
    }
}

#[test]
fn affine_correlated_likelihood_matches_one_existing_batch_update() {
    with_cx(&CancelGate::new(), |cx| {
        let prior = Belief::new(vec![0.5, -1.0], vec![vec![2.0, 0.5], vec![0.5, 1.0]], cx).unwrap();
        let covariance = vec![vec![0.5, 0.3], vec![0.3, 2.0]];
        let measurements = inputs(vec![2.0, -1.0], covariance.clone(), cx);
        let records = vec![
            Observation::new(vec![1.0, 2.0], 1.75, 0.5, "sensor-0").unwrap(),
            Observation::new(vec![-2.0, 1.0], -0.5, 2.0, "sensor-1").unwrap(),
        ]
        .into_iter()
        .map(RobustObservation::available)
        .collect();
        let batch = ObservationBatch::new(records, covariance.clone(), cx).unwrap();
        let oracle = assimilate_observation_batch(&prior, &batch, cx).unwrap();
        let model = |x: &[f64], _: &Cx<'_>| {
            Ok(ModelEvaluation::Value {
                predictions: vec![x[0] + 2.0 * x[1] + 0.25, -2.0 * x[0] + x[1] - 0.5],
                jacobian: vec![vec![1.0, 2.0], vec![-2.0, 1.0]],
            })
        };
        for iterations in [2, 64] {
            let mut policy = settings();
            policy.max_iterations = iterations;
            let result =
                assimilate_correlated_iterated(&prior, &measurements, policy, model, cx).unwrap();
            assert_eq!(result.outcome(), IteratedOutcome::Stationary);
            assert_eq!(result.trace().len(), 2);
            assert_eq!(result.model_evaluations(), 2);
            belief_close(
                result.posterior().unwrap(),
                oracle.posterior().unwrap(),
                1e-11,
            );
        }
        // A constant reference channel still reveals correlated sensor noise.
        // Conditioning epsilon_0 on epsilon_1=2 gives a reading -1 with R=3/4.
        let prior = Belief::scalar(0.0, 1.0).unwrap();
        let reference = inputs(vec![0.0, 2.0], vec![vec![1.0, 0.5], vec![0.5, 1.0]], cx);
        let model = |x: &[f64], _: &Cx<'_>| {
            Ok(ModelEvaluation::Value {
                predictions: vec![x[0], 0.0],
                jacobian: vec![vec![1.0], vec![0.0]],
            })
        };
        let result =
            assimilate_correlated_iterated(&prior, &reference, settings(), model, cx).unwrap();
        assert_eq!(result.outcome(), IteratedOutcome::Stationary);
        close(result.state()[0], -4.0 / 7.0, 1e-12);
        close(
            result.posterior().unwrap().variance(0).unwrap(),
            3.0 / 7.0,
            1e-12,
        );
    });
}

#[test]
fn shared_noise_does_not_average_away_like_independent_sensor_error() {
    with_cx(&CancelGate::new(), |cx| {
        let prior = Belief::scalar(0.0, 100.0).unwrap();
        // Four sensors: independent variance 1 plus common calibration variance 9.
        let covariance = (0..4)
            .map(|i| (0..4).map(|j| if i == j { 10.0 } else { 9.0 }).collect())
            .collect();
        let correlated = inputs(vec![2.0; 4], covariance, cx);
        let diagonal = IndependentMeasurements::new(
            vec![2.0; 4],
            vec![10.0; 4],
            instruments(4),
            "calibration-v1",
        )
        .unwrap();
        let model = |x: &[f64], _: &Cx<'_>| {
            Ok(ModelEvaluation::Value {
                predictions: vec![x[0]; 4],
                jacobian: vec![vec![1.0]; 4],
            })
        };
        let result =
            assimilate_correlated_iterated(&prior, &correlated, settings(), model, cx).unwrap();
        let naive = assimilate_iterated(&prior, &diagonal, settings(), model, cx).unwrap();
        assert_eq!(result.outcome(), IteratedOutcome::Stationary);
        let variance = result.posterior().unwrap().variance(0).unwrap();
        // The average has variance 9 + 1/4; these scalar formulas do not whiten R.
        close(variance, 1.0 / (0.01 + 1.0 / 9.25), 1e-10);
        close(result.state()[0], 200.0 / 109.25, 1e-11);
        let infinite_sensor_floor = 1.0 / (0.01 + 1.0 / 9.0);
        assert!(variance > infinite_sensor_floor);
        assert!(naive.posterior().unwrap().variance(0).unwrap() < infinite_sensor_floor);
    });
}

#[test]
fn diagonal_covariance_agrees_with_the_independent_nonlinear_lane() {
    with_cx(&CancelGate::new(), |cx| {
        let prior = Belief::scalar(1.0, 1.0).unwrap();
        let correlated = inputs(vec![4.0, 1.8], vec![vec![0.25, 0.0], vec![0.0, 0.5]], cx);
        let independent = IndependentMeasurements::new(
            vec![4.0, 1.8],
            vec![0.25, 0.5],
            instruments(2),
            "calibration-v1",
        )
        .unwrap();
        let actual =
            assimilate_correlated_iterated(&prior, &correlated, settings(), nonlinear, cx).unwrap();
        let expected =
            assimilate_iterated(&prior, &independent, settings(), nonlinear, cx).unwrap();
        assert_eq!(actual.outcome(), IteratedOutcome::Stationary);
        assert_eq!(expected.outcome(), IteratedOutcome::Stationary);
        belief_close(
            actual.posterior().unwrap(),
            expected.posterior().unwrap(),
            2e-9,
        );
    });
}

#[test]
fn consistent_sensor_and_covariance_permutation_preserves_the_physical_answer() {
    with_cx(&CancelGate::new(), |cx| {
        let prior = Belief::scalar(1.0, 1.0).unwrap();
        let original = inputs(vec![4.0, 1.8], noise(), cx);
        let mut ids = instruments(2);
        ids.reverse();
        let permuted = CorrelatedMeasurements::new(
            vec![1.8, 4.0],
            vec![vec![0.5, 0.125], vec![0.125, 0.25]],
            ids,
            "calibration-v1",
            cx,
        )
        .unwrap();
        let swapped = |x: &[f64], _: &Cx<'_>| {
            Ok(ModelEvaluation::Value {
                predictions: vec![x[0], x[0] * x[0]],
                jacobian: vec![vec![1.0], vec![2.0 * x[0]]],
            })
        };
        let a =
            assimilate_correlated_iterated(&prior, &original, settings(), nonlinear, cx).unwrap();
        let b = assimilate_correlated_iterated(&prior, &permuted, settings(), swapped, cx).unwrap();
        assert_eq!(a.outcome(), IteratedOutcome::Stationary);
        assert_eq!(b.outcome(), IteratedOutcome::Stationary);
        belief_close(a.posterior().unwrap(), b.posterior().unwrap(), 2e-9);
    });
}

#[test]
fn nonlinear_optimum_matches_an_independently_inverted_covariance_objective() {
    // R^-1 = [[32,-8],[-8,16]]/7. Differentiate the full scalar objective,
    // including its cross-residual term, and bisect its unique root in [1.5,2].
    let gradient = |x: f64| {
        let (r0, r1) = (x * x - 4.0, x - 1.8);
        x - 1.0 + (2.0 * x * (32.0 * r0 - 8.0 * r1) - 8.0 * r0 + 16.0 * r1) / 7.0
    };
    let (mut lo, mut hi) = (1.5, 2.0);
    for _ in 0..80 {
        let midpoint = (lo + hi) * 0.5;
        if gradient(midpoint) > 0.0 {
            hi = midpoint;
        } else {
            lo = midpoint;
        }
    }
    let optimum = 0.5 * (lo + hi);
    with_cx(&CancelGate::new(), |cx| {
        let prior = Belief::scalar(1.0, 1.0).unwrap();
        let measurements = inputs(vec![4.0, 1.8], noise(), cx);
        let result =
            assimilate_correlated_iterated(&prior, &measurements, settings(), nonlinear, cx)
                .unwrap();
        assert_eq!(result.outcome(), IteratedOutcome::Stationary);
        let x = result.state()[0];
        close(x, optimum, 2e-9);
        assert!(gradient(x).abs() <= 1.1e-7);
        let (r0, r1) = (x * x - 4.0, x - 1.8);
        let objective =
            0.5 * (x - 1.0).powi(2) + (16.0 * r0 * r0 - 8.0 * r0 * r1 + 8.0 * r1 * r1) / 7.0;
        close(result.trace().last().unwrap().objective, objective, 1e-12);
        let local_variance = 1.0 / (1.0 + (128.0 * x * x - 32.0 * x + 16.0) / 7.0);
        close(
            result.posterior().unwrap().variance(0).unwrap(),
            local_variance,
            1e-11,
        );
    });
}

#[test]
fn malformed_nonsymmetric_and_nonpositive_covariances_are_refused() {
    with_cx(&CancelGate::new(), |cx| {
        for covariance in [
            vec![],
            vec![vec![1.0, 0.0], vec![0.0]],
            vec![vec![1.0, 1.0], vec![1.0, 1.0]],
            vec![vec![1.0, 2.0], vec![2.0, 1.0]],
            vec![vec![1.0, 0.5], vec![0.25, 1.0]],
            vec![vec![1.0, f64::NAN], vec![f64::NAN, 1.0]],
            vec![vec![0.0, 0.0], vec![0.0, 1.0]],
        ] {
            assert!(
                CorrelatedMeasurements::new(
                    vec![4.0, 1.8],
                    covariance,
                    instruments(2),
                    "calibration-v1",
                    cx,
                )
                .is_err()
            );
        }
        assert!(
            CorrelatedMeasurements::new(
                vec![4.0, 1.8],
                noise(),
                instruments(1),
                "calibration-v1",
                cx,
            )
            .is_err()
        );
    });
}

#[test]
fn cancellation_and_numerical_caps_never_publish_partial_correlated_posteriors() {
    for before_callback in [true, false] {
        let gate = CancelGate::new();
        with_cx(&gate, |cx| {
            let prior = Belief::scalar(1.0, 1.0).unwrap();
            let measurements = inputs(vec![4.0, 1.8], noise(), cx);
            let calls = Cell::new(0);
            if before_callback {
                gate.request();
            }
            let model = |x: &[f64], cx: &Cx<'_>| {
                calls.set(calls.get() + 1);
                gate.request();
                nonlinear(x, cx)
            };
            let result =
                assimilate_correlated_iterated(&prior, &measurements, settings(), model, cx);
            assert!(matches!(result, Err(AssimError::Cancelled { .. })));
            assert_eq!(calls.get(), usize::from(!before_callback));
        });
    }
    with_cx(&CancelGate::new(), |cx| {
        let prior = Belief::scalar(1.0, 1.0).unwrap();
        let measurements = inputs(vec![4.0, 1.8], noise(), cx);
        for (iterations, evaluations, expected) in [
            (0, 256, IteratedOutcome::IterationLimit),
            (64, 1, IteratedOutcome::EvaluationLimit),
        ] {
            let mut policy = settings();
            policy.max_iterations = iterations;
            policy.max_model_evaluations = evaluations;
            let result =
                assimilate_correlated_iterated(&prior, &measurements, policy, nonlinear, cx)
                    .unwrap();
            assert_eq!(result.outcome(), expected);
            assert!(result.posterior().is_none());
            assert_eq!(result.state(), prior.mean());
            assert_eq!(result.model_evaluations(), 1);
        }
    });
}

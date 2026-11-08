//! Robust nonlinear fits against independently solved contaminated-data optima.

use std::cell::Cell;

use fs_assimilate::nonlinear::iterated::{
    HuberLoss, IndependentMeasurements, IteratedOutcome, IteratedSettings, ModelEvaluation,
    assimilate_iterated, fit_huber_iterated,
};
use fs_assimilate::{AssimError, Belief};
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
                kernel_id: 15,
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

fn measurements(readings: Vec<f64>, model_id: &str) -> IndependentMeasurements {
    let count = readings.len();
    IndependentMeasurements::new(
        readings,
        vec![1.0; count],
        (0..count).map(|i| format!("sensor-{i}")).collect(),
        model_id,
    )
    .unwrap()
}

fn location(state: &[f64], _: &Cx<'_>) -> Result<ModelEvaluation, AssimError> {
    Ok(ModelEvaluation::Value {
        predictions: vec![state[0]; 5],
        jacobian: vec![vec![1.0]; 5],
    })
}

fn square(state: &[f64], _: &Cx<'_>) -> Result<ModelEvaluation, AssimError> {
    Ok(ModelEvaluation::Value {
        predictions: vec![state[0] * state[0]; 4],
        jacobian: vec![vec![2.0 * state[0]]; 4],
    })
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn one_contaminated_reading_has_bounded_influence_on_the_location_fit() {
    // With prior N(0,100), four readings at 2 and one at 100, delta=1,
    // the Huber score equation in the optimum's active region is
    // x/100 + 4*(x-2) - 1 = 0. The Gaussian optimum instead is 108/5.01.
    let optimum = 9.0 / 4.01;
    with_cx(&CancelGate::new(), |cx| {
        let prior = Belief::scalar(0.0, 100.0).unwrap();
        let inputs = measurements(vec![2.0, 2.0, 2.0, 2.0, 100.0], "location-v1");
        let loss = HuberLoss::new(1.0).unwrap();
        let robust = fit_huber_iterated(&prior, &inputs, loss, settings(), location, cx).unwrap();
        assert_eq!(robust.outcome(), IteratedOutcome::Stationary);
        close(robust.state()[0], optimum, 1e-8);
        assert_eq!(robust.huber_loss(), Some(loss));
        assert_eq!(robust.downweighted_rows(), &[4]);
        assert!(robust.posterior().is_none());

        let gaussian = assimilate_iterated(&prior, &inputs, settings(), location, cx).unwrap();
        assert_eq!(gaussian.outcome(), IteratedOutcome::Stationary);
        close(gaussian.state()[0], 108.0 / 5.01, 1e-10);
        assert!((robust.state()[0] - 2.0).abs() < 0.25);
        assert!((gaussian.state()[0] - 2.0).abs() > 19.0);
        assert_eq!(gaussian.huber_loss(), None);
        assert!(gaussian.downweighted_rows().is_empty());
        assert!(gaussian.posterior().is_some());
        assert_ne!(robust.receipt_id(), gaussian.receipt_id());
    });
}

#[test]
fn nonlinear_fit_descends_the_actual_loss_to_an_independent_optimum() {
    // For prior N(1,1), h(x)=x², readings [4,4,4,100] and delta=1,
    // the first three residuals are quadratic at the positive optimum and the
    // fourth is in the positive tail: J'(x) = 6*x³ - 25*x - 1.
    // This derivative is strictly increasing and changes sign on [2,2.2].
    let (mut lo, mut hi) = (2.0, 2.2);
    for _ in 0..80 {
        let x = 0.5 * (lo + hi);
        if 6.0 * x * x * x - 25.0 * x - 1.0 > 0.0 {
            hi = x;
        } else {
            lo = x;
        }
    }
    let optimum = 0.5 * (lo + hi);
    with_cx(&CancelGate::new(), |cx| {
        let prior = Belief::scalar(1.0, 1.0).unwrap();
        let inputs = measurements(vec![4.0, 4.0, 4.0, 100.0], "square-v1");
        let loss = HuberLoss::new(1.0).unwrap();
        let result = fit_huber_iterated(&prior, &inputs, loss, settings(), square, cx).unwrap();
        assert_eq!(result.outcome(), IteratedOutcome::Stationary);
        close(result.state()[0], optimum, 1e-8);
        assert!(result.posterior().is_none());
        assert!(result.trace().len() > 2);
        assert_eq!(result.downweighted_rows(), &[3]);
        assert!(
            result
                .trace()
                .windows(2)
                .all(|w| w[1].objective <= w[0].objective)
        );
        assert!(result.trace().last().unwrap().objective < result.trace()[0].objective);
        assert!(result.trace().last().unwrap().gradient_norm <= settings().gradient_tolerance);
        let x = result.state()[0];
        let clean_residual = 4.0 - x * x;
        assert!(clean_residual.abs() < loss.transition());
        let actual_objective = 0.5 * (x - 1.0) * (x - 1.0)
            + 1.5 * clean_residual * clean_residual
            + (100.0 - x * x - 0.5);
        close(
            result.trace().last().unwrap().objective,
            actual_objective,
            1e-12,
        );
    });
}

#[test]
fn quadratic_region_matches_gaussian_state_but_never_publishes_a_posterior() {
    with_cx(&CancelGate::new(), |cx| {
        let prior = Belief::scalar(0.0, 4.0).unwrap();
        let inputs = measurements(vec![0.5, 1.0, 1.5, 1.0, 1.0], "location-v1");
        let gaussian = assimilate_iterated(&prior, &inputs, settings(), location, cx).unwrap();
        let mut previous_receipt = None;
        for transition in [10.0, 20.0] {
            let loss = HuberLoss::new(transition).unwrap();
            let robust =
                fit_huber_iterated(&prior, &inputs, loss, settings(), location, cx).unwrap();
            assert_eq!(robust.outcome(), IteratedOutcome::Stationary);
            close(robust.state()[0], gaussian.state()[0], 1e-12);
            assert_eq!(robust.huber_loss(), Some(loss));
            assert!(robust.posterior().is_none());
            assert_ne!(robust.receipt_id(), gaussian.receipt_id());
            if let Some(previous) = previous_receipt.replace(robust.receipt_id().to_owned()) {
                assert_ne!(previous, robust.receipt_id());
            }
            let replay =
                fit_huber_iterated(&prior, &inputs, loss, settings(), location, cx).unwrap();
            assert_eq!(robust, replay);
        }

        // Even an initially stationary fit has no Gaussian posterior authority.
        let inputs = measurements(vec![0.0; 5], "location-v1");
        let result = fit_huber_iterated(
            &prior,
            &inputs,
            HuberLoss::new(1.0).unwrap(),
            settings(),
            location,
            cx,
        )
        .unwrap();
        assert_eq!(result.outcome(), IteratedOutcome::Stationary);
        assert_eq!(result.model_evaluations(), 1);
        assert!(result.posterior().is_none());
    });
}

#[test]
fn numerical_limits_report_unresolved_huber_fits_without_posterior() {
    with_cx(&CancelGate::new(), |cx| {
        let prior = Belief::scalar(1.0, 1.0).unwrap();
        let inputs = measurements(vec![4.0, 4.0, 4.0, 100.0], "square-v1");
        let loss = HuberLoss::new(1.0).unwrap();
        for (iterations, evaluations, expected, consumed) in [
            (0, 256, IteratedOutcome::IterationLimit, 1),
            (64, 1, IteratedOutcome::EvaluationLimit, 1),
            (1, 256, IteratedOutcome::IterationLimit, 2),
        ] {
            let mut policy = settings();
            policy.max_iterations = iterations;
            policy.max_model_evaluations = evaluations;
            let result = fit_huber_iterated(&prior, &inputs, loss, policy, square, cx).unwrap();
            assert_eq!(result.outcome(), expected);
            assert_eq!(result.model_evaluations(), consumed);
            assert_eq!(result.huber_loss(), Some(loss));
            assert!(result.posterior().is_none());
        }
    });
}

#[test]
fn cancellation_before_or_during_model_evaluation_does_not_publish_a_fit() {
    for pre_cancelled in [true, false] {
        let gate = CancelGate::new();
        let calls = Cell::new(0);
        if pre_cancelled {
            gate.request();
        }
        with_cx(&gate, |cx| {
            let prior = Belief::scalar(1.0, 1.0).unwrap();
            let original = prior.clone();
            let inputs = measurements(vec![4.0, 4.0, 4.0, 100.0], "square-v1");
            let result = fit_huber_iterated(
                &prior,
                &inputs,
                HuberLoss::new(1.0).unwrap(),
                settings(),
                |state, cx| {
                    calls.set(calls.get() + 1);
                    gate.request();
                    square(state, cx)
                },
                cx,
            );
            assert!(matches!(result, Err(AssimError::Cancelled { .. })));
            assert_eq!(calls.get(), usize::from(!pre_cancelled));
            assert_eq!(prior, original);
        });
    }
}

//! G0/G3/G4 nonlinear finite-difference linearization and demotion tests.

use std::cell::Cell;

use fs_assimilate::{
    AssimError, Belief, FiniteDifferenceSettings, LinearizationDisposition,
    NonlinearObservationSpec, Observation, assimilate, assimilate_nonlinear_fd,
    linearize_nonlinear_fd,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};

const TEST_STREAM: StreamKey = StreamKey {
    seed: 0x4E4F_4E4C_494E_4541,
    kernel_id: 0x4644,
    tile: 0,
    iteration: 0,
};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            TEST_STREAM,
            Budget::INFINITE,
            ExecMode::Deterministic,
        )
        .with_time_source(&clock);
        f(&cx)
    });
    assert!(pool.stats().quiescent());
    result
}

fn settings(relative_step: f64, gate: f64) -> FiniteDifferenceSettings {
    FiniteDifferenceSettings::new(relative_step, gate).expect("valid settings")
}

fn spec(
    reading: f64,
    noise_var: f64,
    model_id: &str,
    settings: FiniteDifferenceSettings,
) -> NonlinearObservationSpec {
    NonlinearObservationSpec::new(reading, noise_var, "nonlinear-sensor-a", model_id, settings)
        .expect("valid nonlinear observation")
}

#[test]
fn linear_model_reduces_exactly_to_the_existing_scalar_updater() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let prior = Belief::scalar(2.0, 1.5).expect("valid prior");
        let nonlinear_spec = spec(5.5, 0.5, "affine-two-x-plus-one-v1", settings(0.25, 10.0));
        let evaluations = Cell::new(0_u32);
        let result = assimilate_nonlinear_fd(
            &prior,
            &nonlinear_spec,
            |state| {
                evaluations.set(evaluations.get() + 1);
                2.0 * state[0] + 1.0
            },
            cx,
        )
        .expect("admitted nonlinear update");

        assert_eq!(evaluations.get(), 3);
        assert_eq!(
            result.linearization().disposition(),
            LinearizationDisposition::Admitted
        );
        assert_eq!(result.linearization().prediction(), 5.0);
        assert_eq!(result.linearization().innovation(), 0.5);
        assert_eq!(result.linearization().jacobian(), &[2.0]);

        let expected_observation = Observation::new(vec![2.0], 4.5, 0.5, "nonlinear-sensor-a")
            .expect("valid affine observation");
        assert_eq!(
            result.linearization().observation(),
            Some(&expected_observation)
        );
        let expected = assimilate(&prior, &expected_observation, cx).expect("direct linear update");
        assert_eq!(result.posterior(), Some(&expected));
    });
}

#[test]
fn probes_retain_actual_states_and_receipts_bind_model_identity() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let prior = Belief::diagonal(vec![2.0, -1.0], &[0.5, 2.0], cx).expect("valid prior");
        let first_spec = spec(1.25, 0.25, "quadratic-plus-linear-v1", settings(0.25, 10.0));
        let predict = |state: &[f64]| state[0] * state[0] + 3.0 * state[1];
        let first =
            linearize_nonlinear_fd(&prior, &first_spec, predict, cx).expect("linearization");
        let replay =
            linearize_nonlinear_fd(&prior, &first_spec, predict, cx).expect("deterministic replay");

        assert_eq!(first, replay);
        assert_eq!(first.prediction(), 1.0);
        assert_eq!(first.jacobian(), &[4.0, 3.0]);
        assert_eq!(first.probes().len(), 2);
        assert_eq!(first.probes()[0].component(), 0);
        assert_eq!(first.probes()[0].nominal_step(), 0.5);
        assert_eq!(first.probes()[0].plus_state(), 2.5);
        assert_eq!(first.probes()[0].minus_state(), 1.5);
        assert_eq!(first.probes()[0].plus_prediction(), 3.25);
        assert_eq!(first.probes()[0].minus_prediction(), -0.75);
        assert_eq!(first.probes()[0].denominator(), 1.0);
        assert_eq!(first.probes()[1].component(), 1);
        assert_eq!(first.probes()[1].nominal_step(), 0.25);
        assert_eq!(first.probes()[1].plus_state(), -0.75);
        assert_eq!(first.probes()[1].minus_state(), -1.25);
        assert_eq!(first.probes()[1].plus_prediction(), 1.75);
        assert_eq!(first.probes()[1].minus_prediction(), 0.25);
        assert!(
            first
                .receipt_id()
                .starts_with("nonlinear-fd-linearization:v1:")
        );

        let renamed_spec = spec(1.25, 0.25, "quadratic-plus-linear-v2", settings(0.25, 10.0));
        let renamed =
            linearize_nonlinear_fd(&prior, &renamed_spec, predict, cx).expect("renamed model");
        assert_eq!(first.jacobian(), renamed.jacobian());
        assert_ne!(first.receipt_id(), renamed.receipt_id());
    });
}

#[test]
fn large_innovation_is_demoted_without_observation_or_posterior() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let prior = Belief::scalar(0.0, 1.0).expect("valid prior");
        let nonlinear_spec = spec(100.0, 1.0, "square-v1", settings(0.25, 3.0));
        let result =
            assimilate_nonlinear_fd(&prior, &nonlinear_spec, |state| state[0] * state[0], cx)
                .expect("audited demotion");

        assert_eq!(
            result.linearization().disposition(),
            LinearizationDisposition::DemotedLargeInnovation
        );
        assert_eq!(result.linearization().standardized_innovation(), 100.0);
        assert!(result.linearization().observation().is_none());
        assert!(result.posterior().is_none());
    });
}

#[test]
fn checked_constructors_refuse_invalid_policy_and_identity() {
    assert_eq!(
        FiniteDifferenceSettings::new(0.0, 3.0),
        Err(AssimError::InvalidLinearizationParameter {
            parameter: "relative finite-difference step"
        })
    );
    assert_eq!(
        FiniteDifferenceSettings::new(0.25, f64::INFINITY),
        Err(AssimError::InvalidLinearizationParameter {
            parameter: "maximum standardized innovation"
        })
    );
    assert_eq!(
        NonlinearObservationSpec::new(f64::NAN, 1.0, "sensor", "model", settings(0.25, 3.0)),
        Err(AssimError::NonFiniteObservationValue)
    );
    assert!(matches!(
        NonlinearObservationSpec::new(0.0, 1.0, "sensor", "unknown", settings(0.25, 3.0)),
        Err(AssimError::InvalidIdentity {
            field: "nonlinear_model",
            ..
        })
    ));
}

#[test]
fn nonfinite_model_and_jacobian_paths_are_typed_refusals() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let prior = Belief::scalar(0.0, 1.0).expect("valid prior");
        let nonlinear_spec = spec(0.0, 1.0, "nonfinite-probe-v1", settings(0.25, 3.0));
        assert_eq!(
            linearize_nonlinear_fd(&prior, &nonlinear_spec, |_| f64::NAN, cx),
            Err(AssimError::NonFiniteModelPrediction {
                evaluation: "base",
                component: None
            })
        );

        assert_eq!(
            linearize_nonlinear_fd(
                &prior,
                &nonlinear_spec,
                |state| if state[0] > 0.0 { f64::NAN } else { 0.0 },
                cx
            ),
            Err(AssimError::NonFiniteModelPrediction {
                evaluation: "plus",
                component: Some(0)
            })
        );

        assert_eq!(
            linearize_nonlinear_fd(
                &prior,
                &nonlinear_spec,
                |state| {
                    if state[0] > 0.0 {
                        f64::MAX
                    } else if state[0] < 0.0 {
                        -f64::MAX
                    } else {
                        0.0
                    }
                },
                cx
            ),
            Err(AssimError::NonFiniteObservationJacobian { component: 0 })
        );
    });
}

#[test]
fn unrepresentable_step_is_refused_before_probe_publication() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let prior = Belief::scalar(f64::MAX, 1.0).expect("finite prior");
        let nonlinear_spec = spec(0.0, 1.0, "overflow-step-v1", settings(1.0, 3.0));
        assert_eq!(
            linearize_nonlinear_fd(&prior, &nonlinear_spec, |state| state[0], cx),
            Err(AssimError::FiniteDifferenceStepUnrepresentable { component: 0 })
        );

        let prior = Belief::scalar(0.0, 1.0).expect("finite prior");
        let nonlinear_spec = spec(0.0, 1.0, "overflow-denominator-v1", settings(f64::MAX, 3.0));
        assert_eq!(
            linearize_nonlinear_fd(&prior, &nonlinear_spec, |state| state[0], cx),
            Err(AssimError::FiniteDifferenceStepUnrepresentable { component: 0 })
        );
    });
}

#[test]
fn precancelled_linearization_publishes_nothing() {
    let gate = CancelGate::new();
    gate.request();
    with_cx(&gate, |cx| {
        let prior = Belief::scalar(0.0, 1.0).expect("valid prior");
        let nonlinear_spec = spec(0.0, 1.0, "cancelled-model-v1", settings(0.25, 3.0));
        let error = linearize_nonlinear_fd(&prior, &nonlinear_spec, |state| state[0], cx)
            .expect_err("pre-cancelled call must refuse");
        assert!(matches!(
            error,
            AssimError::Cancelled {
                phase: "nonlinear preflight",
                completed: 0,
                ..
            }
        ));
    });
}

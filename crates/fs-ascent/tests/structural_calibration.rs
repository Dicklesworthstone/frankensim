//! G3/G4/G5: physical calibration, total gradients, bounded refusal and replay.
#[allow(dead_code)]
#[path = "../examples/structural_calibration.rs"]
mod structural;

use fs_ascent::transient::structural::{
    StructuralTransientError, StructuralTransientFamily, StructuralTransientModel,
    StructuralTransientStudy, evaluate_structural_transient,
};
use fs_ascent::{SqpError, SqpStop};
use fs_solver::FlexiblePreconditioner;
use fs_time::galpha::SecondOrderProblem;
use fs_time::galpha::second_order_adjoint::trajectory::{
    StructuralRecordingStatus, StructuralTrajectoryError, StructuralTrajectoryModel,
};
use fs_time::stiff::IdentityPreconditioner;
use std::cell::Cell;
use structural::{START, STEPS, StructuralReadings, TRUTH, config};

fn check_gradient(data: &StructuralReadings, point: &[f64; 5], actual: &[f64]) {
    for j in 0..5 {
        let (mut plus, mut minus) = (*point, *point);
        let h = 2e-5;
        plus[j] += h;
        minus[j] -= h;
        let expected = (data.dense_loss(&plus) - data.dense_loss(&minus)) / (2.0 * h);
        assert!(
            (actual[j] - expected).abs() < 2e-7 * (1.0 + expected.abs()),
            "parameter {j}: adjoint={} independent dense FD={expected}",
            actual[j]
        );
    }
}

#[test]
fn physical_initialization_and_all_parameter_routes_match_independent_dense_gradients() {
    let data = StructuralReadings::synthetic();
    let model = data.instantiate(&TRUTH).unwrap();
    let initial = model.initial_state();
    let (mut inertia, mut damping, mut spring, mut force) =
        ([0.0; 2], [0.0; 2], [0.0; 2], [0.0; 2]);
    model.mass_apply(&initial.a, &mut inertia);
    model.damping_apply(&initial.v, &mut damping);
    model.internal_force(&initial.q, &mut spring);
    model.forcing(0.0, &mut force).unwrap();
    for i in 0..2 {
        assert!((inertia[i] + damping[i] + spring[i] - force[i]).abs() < 1e-12);
    }
    // Internal coupling cancels in net force; damping dissipates physical power.
    assert!((spring.iter().sum::<f64>() - 0.8 * initial.q[0] - 1.1 * initial.q[1]).abs() < 1e-12);
    let dissipation = damping[0] * initial.v[0] + damping[1] * initial.v[1];
    let expected = 0.08 * initial.v[0].powi(2)
        + 0.12 * initial.v[1].powi(2)
        + TRUTH[1] * (initial.v[0] - initial.v[1]).powi(2);
    assert!(dissipation > 0.0 && (dissipation - expected).abs() < 1e-12);
    let mut later = [0.0; 2];
    model.forcing(0.73, &mut later).unwrap();
    assert_ne!(force, later);

    for point in [START, [1.6, 0.12, 1.4, 0.35, 0.1]] {
        let got = evaluate_structural_transient(
            &data,
            &config(),
            &point,
            &IdentityPreconditioner,
            &mut || false,
        )
        .unwrap()
        .unwrap();
        assert!((got.value - data.dense_loss(&point)).abs() < 1e-9);
        check_gradient(&data, &point, &got.gradient);
        assert_eq!(got.observations, data.indices.len());
        assert_eq!(got.forward.status, StructuralRecordingStatus::ReachedEnd);
        assert_eq!(got.forward.advanced, STEPS);
        assert_eq!(got.final_state.steps, STEPS);
        assert_eq!(got.final_state.history, Vec::new());
        assert!(got.replayed_steps >= STEPS);
        assert!(got.peak_checkpoints <= config().replay.checkpoints);
        let repeated = evaluate_structural_transient(
            &data,
            &config(),
            &point,
            &IdentityPreconditioner,
            &mut || false,
        )
        .unwrap()
        .unwrap();
        assert_eq!(got, repeated);
    }
}

#[test]
fn sqp_recovers_five_parameters_with_bit_identical_split_and_clone_continuation() {
    let data = StructuralReadings::synthetic();
    let adjoint = IdentityPreconditioner;
    let mut full =
        StructuralTransientStudy::new(&data, &START, config(), &adjoint, &mut || false).unwrap();
    let mut split = full.clone();
    let initial = full.accepted().value;
    let report = full.run(1e-7, 100, 1000, &mut || false).unwrap();
    assert_eq!(report.stop, SqpStop::Converged, "{report:?}");
    for (actual, expected) in full.optimizer().point().iter().zip(TRUTH) {
        assert!((actual - expected).abs() < 2e-5, "{actual} != {expected}");
    }
    assert!(full.accepted().value < initial * 1e-10);
    assert_eq!(
        split.run(1e-7, 3, 1000, &mut || false).unwrap().stop,
        SqpStop::IterationLimit
    );
    let mut fork = split.clone();
    for _ in 0..100 {
        if split.run(1e-7, 1, 1000, &mut || false).unwrap().stop != SqpStop::IterationLimit {
            break;
        }
    }
    assert_eq!(
        fork.run(1e-7, 100, 1000, &mut || false).unwrap().stop,
        SqpStop::Converged
    );
    for other in [&split, &fork] {
        assert_eq!(full.optimizer().point(), other.optimizer().point());
        assert_eq!(full.optimizer().history(), other.optimizer().history());
        assert_eq!(
            full.optimizer().evaluations(),
            other.optimizer().evaluations()
        );
        assert_eq!(full.accepted(), other.accepted());
    }
    assert_eq!(full.accepted().point, full.optimizer().point());
    assert_eq!(full.accepted().value, full.optimizer().sample().f);
    assert_eq!(full.accepted().gradient, full.optimizer().sample().gradient);
}

#[derive(Default)]
struct ProbePreconditioner {
    calls: Cell<usize>,
    fail: Cell<bool>,
}
impl FlexiblePreconditioner for ProbePreconditioner {
    fn apply(&self, _: usize, residual: &[f64], output: &mut [f64]) {
        self.calls.set(self.calls.get() + 1);
        if self.fail.get() {
            output.fill(f64::NAN);
        } else {
            output.copy_from_slice(residual);
        }
    }
}

#[test]
fn cancellation_and_adjoint_refusal_preserve_accepted_physics_and_charge_attempts() {
    let data = StructuralReadings::synthetic();
    let adjoint = ProbePreconditioner::default();
    let mut study =
        StructuralTransientStudy::new(&data, &START, config(), &adjoint, &mut || false).unwrap();
    let before = study.accepted().clone();
    let calls = adjoint.calls.get();
    assert!(
        calls > 0,
        "the supplied adjoint preconditioner must be used"
    );
    assert_eq!(
        study.run(1e-7, 10, 1, &mut || false).unwrap().stop,
        SqpStop::EvaluationLimit
    );
    assert_eq!(adjoint.calls.get(), calls);
    assert_eq!(study.accepted(), &before);
    assert!(matches!(
        study.run(1e-7, 1, 1000, &mut || adjoint.calls.get() > calls),
        Err(SqpError::Cancelled)
    ));
    assert!(study.optimizer().evaluations() > 1);
    assert_eq!(study.accepted(), &before);
    assert_eq!(study.optimizer().point(), before.point);
    let spent = study.optimizer().evaluations();
    adjoint.fail.set(true);
    assert!(matches!(
        study.run(1e-7, 1, 1000, &mut || false),
        Err(SqpError::Evaluation(StructuralTransientError::Trajectory(
            _
        )))
    ));
    assert!(study.optimizer().evaluations() > spent);
    assert_eq!(study.accepted(), &before);
    adjoint.fail.set(false);
    assert_eq!(
        study.run(1e-7, 100, 1000, &mut || false).unwrap().stop,
        SqpStop::Converged
    );
    assert!(study.accepted().value < before.value * 1e-10);
}

struct CountedFamily {
    data: StructuralReadings,
    instances: Cell<usize>,
}
impl StructuralTransientFamily for CountedFamily {
    type Model = structural::StructuralModel;
    fn bounds(&self) -> &[[f64; 2]] {
        self.data.bounds()
    }
    fn sample_indices(&self) -> &[usize] {
        self.data.sample_indices()
    }
    fn instantiate(&self, point: &[f64]) -> Result<Self::Model, String> {
        self.instances.set(self.instances.get() + 1);
        self.data.instantiate(point)
    }
}

#[test]
fn admission_and_forward_record_checkpoint_replay_limits_refuse_incomplete_evaluations() {
    let mut data = CountedFamily {
        data: StructuralReadings::synthetic(),
        instances: Cell::new(0),
    };
    let evaluate = |data: &CountedFamily, cfg: &_, point: &[f64]| {
        evaluate_structural_transient(data, cfg, point, &IdentityPreconditioner, &mut || false)
    };
    assert!(
        evaluate(&data, &config(), &[3.0, 0.2, 1.0, 0.4, 0.0])
            .unwrap()
            .is_none()
    );
    assert!(evaluate(&data, &config(), &[f64::NAN, 0.2, 1.0, 0.4, 0.0]).is_err());
    let mut cfg = config();
    cfg.max_samples = 1;
    assert!(matches!(
        evaluate(&data, &cfg, &START),
        Err(StructuralTransientError::Invalid(_))
    ));
    cfg = config();
    cfg.max_kkt_dimension = 14;
    assert!(matches!(
        evaluate(&data, &cfg, &START),
        Err(StructuralTransientError::Invalid(_))
    ));
    cfg = config();
    cfg.recording.adjoint.restart = 0;
    assert!(matches!(
        evaluate(&data, &cfg, &START),
        Err(StructuralTransientError::Invalid(_))
    ));
    data.data.indices.swap(0, 6);
    assert!(matches!(
        evaluate(&data, &config(), &START),
        Err(StructuralTransientError::Invalid(_))
    ));
    data.data.indices.swap(0, 6);
    assert_eq!(
        data.instances.get(),
        0,
        "admission precedes model factory work"
    );

    cfg = config();
    cfg.max_forward_steps = 0;
    assert!(matches!(
        evaluate(&data, &cfg, &START),
        Err(StructuralTransientError::ForwardStopped(
            StructuralRecordingStatus::StepLimit
        ))
    ));
    cfg = config();
    cfg.max_records = 0;
    assert!(matches!(
        evaluate(&data, &cfg, &START),
        Err(StructuralTransientError::ForwardStopped(
            StructuralRecordingStatus::RecordLimit
        ))
    ));
    cfg = config();
    cfg.replay.checkpoints = 0;
    assert!(matches!(
        evaluate(&data, &cfg, &START),
        Err(StructuralTransientError::Trajectory(
            StructuralTrajectoryError::CheckpointLimit { .. }
        ))
    ));
    cfg = config();
    cfg.replay.forward_steps = 0;
    assert!(matches!(
        evaluate(&data, &cfg, &START),
        Err(StructuralTransientError::Trajectory(
            StructuralTrajectoryError::ReplayLimit
        ))
    ));
}

#[test]
fn zero_step_observations_include_consistent_acceleration_and_direct_sensor_derivatives() {
    let mut data = StructuralReadings::synthetic();
    data.indices = vec![0; 6];
    let mut cfg = config();
    cfg.recording.steps = 0;
    cfg.max_forward_steps = 0;
    cfg.max_records = 0;
    cfg.replay.checkpoints = 0;
    cfg.replay.forward_steps = 0;
    let adjoint = ProbePreconditioner::default();
    let got = evaluate_structural_transient(&data, &cfg, &START, &adjoint, &mut || false)
        .unwrap()
        .unwrap();
    assert!((got.value - data.dense_loss(&START)).abs() < 1e-13);
    check_gradient(&data, &START, &got.gradient);
    assert!((got.gradient[4] + 0.4).abs() < 1e-13);
    assert!(got.gradient[..4].iter().all(|g| g.abs() > 1e-3));
    assert_eq!(
        got.final_state,
        data.instantiate(&START).unwrap().initial_state().clone()
    );
    assert_eq!(got.observations, 6);
    assert_eq!((got.replayed_steps, got.peak_checkpoints), (0, 0));
    assert_eq!(adjoint.calls.get(), 0);
}

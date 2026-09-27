//! Full stiff thermal-calibration consumer through the existing SQP engine.
#[allow(dead_code)]
#[path = "../examples/imex_calibration.rs"]
mod thermal;

use fs_ascent::transient::imex::{
    ImexTransientError, ImexTransientFamily, ImexTransientStudy, evaluate_imex_transient,
};
use fs_ascent::{SqpError, SqpStop};
use fs_solver::{FlexiblePreconditioner, LinearOp};
use fs_time::stiff::IdentityPreconditioner;
use fs_time::stiff::adjoint::ImexVjp;
use fs_time::stiff::adjoint::trajectory::{ImexRecordingStatus, ImexTrajectoryError};
use std::cell::Cell;
use thermal::{HeatReadings, START, TRUTH, config};

#[test]
fn nonsymmetric_heat_total_gradient_matches_independent_dense_finite_differences() {
    let data = HeatReadings::synthetic();
    let model = data.instantiate(&TRUTH).unwrap();
    let (u, seed) = ([1.1, -0.3], [0.4, 0.7]);
    let (mut lu, mut lt_seed, mut source) = ([0.0; 2], [0.0; 2], [0.0; 2]);
    model.apply(&u, &mut lu);
    model.apply_transpose(&seed, &mut lt_seed);
    model.nonlinear(&u, &mut source);
    let dot = |a: [f64; 2], b: [f64; 2]| a[0] * b[0] + a[1] * b[1];
    assert!((dot(lu, seed) - dot(u, lt_seed)).abs() < 1e-12);
    // Internal conduction cancels in the physical energy balance C_i*T_i'.
    let energy_rate = 0.01 * (lu[0] + source[0]) + 2.0 * (lu[1] + source[1]);
    assert!((energy_rate - (TRUTH[1] - 0.1 * u[0] - 0.2 * u[1])).abs() < 1e-12);
    model.apply(&[1.0, 0.0], &mut lu);
    model.apply(&[0.0, 1.0], &mut lt_seed);
    assert!((lu[1] - lt_seed[0]).abs() > 100.0);

    let (primal, adjoint) = (IdentityPreconditioner, IdentityPreconditioner);
    for point in [START, [1.6, 2.2, -0.1, 1.1]] {
        let got =
            evaluate_imex_transient(&data, &config(), &point, &primal, &adjoint, &mut || false)
                .unwrap()
                .unwrap();
        assert!((got.value - data.dense_loss(&point)).abs() < 1e-10);
        assert_eq!(got.observations, data.indices.len());
        assert_eq!(got.forward.status, ImexRecordingStatus::ReachedEnd);
        assert_eq!(got.forward.advanced, thermal::STEPS);
        assert!(got.replayed_steps >= thermal::STEPS);
        assert!(got.peak_checkpoints <= config().replay.checkpoints);
        for j in 0..4 {
            let (mut plus, mut minus) = (point, point);
            let h = 1e-5;
            plus[j] += h;
            minus[j] -= h;
            let fd = (data.dense_loss(&plus) - data.dense_loss(&minus)) / (2.0 * h);
            assert!(
                (got.gradient[j] - fd).abs() < 2e-7 * (1.0 + fd.abs()),
                "parameter {j}: adjoint={} dense FD={fd}",
                got.gradient[j]
            );
        }
        let replay =
            evaluate_imex_transient(&data, &config(), &point, &primal, &adjoint, &mut || false)
                .unwrap()
                .unwrap();
        assert_eq!(got, replay);
    }
}

#[test]
fn native_sqp_recovers_four_thermal_parameters_and_continues_bit_identically() {
    let data = HeatReadings::synthetic();
    let (primal, adjoint) = (IdentityPreconditioner, IdentityPreconditioner);
    let mut full =
        ImexTransientStudy::new(&data, &START, config(), &primal, &adjoint, &mut || false).unwrap();
    let mut split =
        ImexTransientStudy::new(&data, &START, config(), &primal, &adjoint, &mut || false).unwrap();
    let initial = full.accepted().value;
    let report = full.run(1e-7, 100, 1000, &mut || false).unwrap();
    assert_eq!(report.stop, SqpStop::Converged, "{report:?}");
    for (actual, expected) in full.optimizer().point().iter().zip(TRUTH) {
        assert!((actual - expected).abs() < 2e-5, "{actual} != {expected}");
    }
    assert!(full.accepted().value < initial * 1e-10);
    for _ in 0..100 {
        if split.run(1e-7, 1, 1000, &mut || false).unwrap().stop != SqpStop::IterationLimit {
            break;
        }
    }
    assert_eq!(full.optimizer().point(), split.optimizer().point());
    assert_eq!(full.optimizer().history(), split.optimizer().history());
    assert_eq!(
        full.optimizer().evaluations(),
        split.optimizer().evaluations()
    );
    assert_eq!(full.accepted(), split.accepted());
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
fn cancelled_or_failed_adjoint_keeps_accepted_physics_and_counts_spent_trials() {
    let data = HeatReadings::synthetic();
    let (primal, adjoint) = (
        ProbePreconditioner::default(),
        ProbePreconditioner::default(),
    );
    let mut study =
        ImexTransientStudy::new(&data, &START, config(), &primal, &adjoint, &mut || false).unwrap();
    let before = study.accepted().clone();
    let calls = (primal.calls.get(), adjoint.calls.get());
    assert!(
        calls.0 > 0 && calls.1 > 0,
        "both preconditioners must be used"
    );
    assert_eq!(
        study.run(1e-7, 10, 1, &mut || false).unwrap().stop,
        SqpStop::EvaluationLimit
    );
    assert_eq!((primal.calls.get(), adjoint.calls.get()), calls);
    assert_eq!(study.accepted(), &before);
    assert!(matches!(
        study.run(1e-7, 1, 100, &mut || adjoint.calls.get() > calls.1),
        Err(SqpError::Cancelled)
    ));
    assert!(study.optimizer().evaluations() > 1);
    assert_eq!(study.accepted(), &before);
    assert_eq!(study.optimizer().point(), before.point);
    let spent = study.optimizer().evaluations();
    adjoint.fail.set(true);
    assert!(matches!(
        study.run(1e-7, 1, 100, &mut || false),
        Err(SqpError::Evaluation(ImexTransientError::Trajectory(
            ImexTrajectoryError::Step(_)
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
    data: HeatReadings,
    instances: Cell<usize>,
}
impl ImexTransientFamily for CountedFamily {
    type Model = thermal::HeatModel;
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
fn invalid_trials_and_exhausted_work_refuse_without_publishing_penalties() {
    let mut data = CountedFamily {
        data: HeatReadings::synthetic(),
        instances: Cell::new(0),
    };
    let (primal, adjoint) = (IdentityPreconditioner, IdentityPreconditioner);
    let evaluate = |data: &CountedFamily, cfg: &_, point: &[f64]| {
        evaluate_imex_transient(data, cfg, point, &primal, &adjoint, &mut || false)
    };
    assert!(
        evaluate(&data, &config(), &[3.0, 1.1, 0.0, 1.0])
            .unwrap()
            .is_none()
    );
    assert!(evaluate(&data, &config(), &[f64::NAN, 1.1, 0.0, 1.0]).is_err());
    let mut cfg = config();
    cfg.max_samples = 1;
    assert!(matches!(
        evaluate(&data, &cfg, &START),
        Err(ImexTransientError::Invalid(_))
    ));
    cfg = config();
    cfg.max_kkt_dimension = 11;
    assert!(matches!(
        evaluate(&data, &cfg, &START),
        Err(ImexTransientError::Invalid(_))
    ));
    cfg = config();
    cfg.solve.restart = 0;
    assert!(matches!(
        evaluate(&data, &cfg, &START),
        Err(ImexTransientError::Invalid(_))
    ));
    data.data.indices.swap(0, 3);
    assert!(matches!(
        evaluate(&data, &config(), &START),
        Err(ImexTransientError::Invalid(_))
    ));
    data.data.indices.swap(0, 3);
    assert_eq!(
        data.instances.get(),
        0,
        "admission must precede the model factory"
    );

    cfg = config();
    cfg.max_forward_steps = 0;
    assert!(matches!(
        evaluate(&data, &cfg, &START),
        Err(ImexTransientError::ForwardStopped(
            ImexRecordingStatus::StepLimit
        ))
    ));
    cfg = config();
    cfg.max_records = 0;
    assert!(matches!(
        evaluate(&data, &cfg, &START),
        Err(ImexTransientError::ForwardStopped(
            ImexRecordingStatus::RecordLimit
        ))
    ));
    cfg = config();
    cfg.replay.checkpoints = 0;
    assert!(matches!(
        evaluate(&data, &cfg, &START),
        Err(ImexTransientError::Trajectory(
            ImexTrajectoryError::CheckpointLimit { .. }
        ))
    ));
    cfg = config();
    cfg.replay.forward_steps = 0;
    assert!(matches!(
        evaluate(&data, &cfg, &START),
        Err(ImexTransientError::Trajectory(
            ImexTrajectoryError::ReplayLimit
        ))
    ));
}

#[test]
fn repeated_initial_sensors_need_no_steps_or_checkpoints() {
    let mut data = HeatReadings::synthetic();
    data.indices = vec![0, 0];
    let mut cfg = config();
    cfg.recording.steps = 0;
    cfg.max_forward_steps = 0;
    cfg.max_records = 0;
    cfg.replay.checkpoints = 0;
    cfg.replay.forward_steps = 0;
    let (primal, adjoint) = (
        ProbePreconditioner::default(),
        ProbePreconditioner::default(),
    );
    let point = [1.1, 1.4, -0.1, 1.6];
    let got = evaluate_imex_transient(&data, &cfg, &point, &primal, &adjoint, &mut || false)
        .unwrap()
        .unwrap();
    assert!((got.value - 0.5 * (0.55_f64.powi(2) + 0.25_f64.powi(2))).abs() < 1e-14);
    for (actual, expected) in got.gradient.iter().zip([0.0, 0.0, 0.3, 0.55]) {
        assert!((actual - expected).abs() < 1e-14);
    }
    assert_eq!(got.final_state, [1.6, 0.0]);
    assert_eq!((got.replayed_steps, got.peak_checkpoints), (0, 0));
    assert_eq!((primal.calls.get(), adjoint.calls.get()), (0, 0));
}

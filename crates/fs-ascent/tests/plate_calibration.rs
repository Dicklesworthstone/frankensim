//! G3/G4/G5: spatial plate calibration, total gradients and accepted-state replay.
#[allow(dead_code)]
#[path = "../examples/plate_calibration.rs"]
mod calibration;

use calibration::{PlateCalibration, START, TRUTH};
use fs_ascent::{SqpError, SqpState, SqpStop};

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn same_accepted(a: &SqpState, b: &SqpState) {
    assert_eq!(bits(a.point()), bits(b.point()));
    assert_eq!(a.sample().f.to_bits(), b.sample().f.to_bits());
    for (left, right) in [
        (&a.sample().gradient, &b.sample().gradient),
        (&a.sample().ce, &b.sample().ce),
        (&a.sample().ci, &b.sample().ci),
        (&a.sample().je, &b.sample().je),
        (&a.sample().ji, &b.sample().ji),
    ] {
        assert_eq!(bits(left), bits(right));
    }
    assert_eq!(bits(a.history()), bits(b.history()));
    assert_eq!(a.iterations(), b.iterations());
}

#[test]
fn complete_plate_gradient_matches_independent_dense_forward_differences() {
    let data = PlateCalibration::synthetic().unwrap();
    let truth = data.evaluate(&TRUTH, &mut || false).unwrap().unwrap();
    assert!(truth.f < 1e-12, "truth objective {}", truth.f);
    assert!(data.dense_loss(&TRUTH).unwrap() < 1e-20);
    let midpoint = std::array::from_fn(|i| 0.5 * (START[i] + TRUTH[i]));
    for point in [START, midpoint] {
        let got = data.evaluate(&point, &mut || false).unwrap().unwrap();
        let reference = data.dense_loss(&point).unwrap();
        assert!((got.f - reference).abs() < 1e-8 * (1.0 + reference.abs()));
        for j in 0..3 {
            let (mut plus, mut minus) = (point, point);
            let h = 2e-5 * (1.0 + point[j].abs());
            plus[j] += h;
            minus[j] -= h;
            let fd =
                (data.dense_loss(&plus).unwrap() - data.dense_loss(&minus).unwrap()) / (2.0 * h);
            assert!(
                (got.gradient[j] - fd).abs() < 2e-6 * (1.0 + fd.abs()),
                "parameter {j}: adjoint={} independent dense FD={fd}",
                got.gradient[j]
            );
        }
    }
}

#[test]
fn spatial_plate_parameters_recover_with_exact_split_and_clone_continuation() {
    let data = PlateCalibration::synthetic().unwrap();
    let mut evaluate = |point: &[f64]| data.evaluate(point, &mut || false);
    let mut full = SqpState::try_new(&START, 9, &mut evaluate, None).unwrap();
    let mut split = full.clone();
    let initial = full.sample().f;
    let report = full.try_run(&mut evaluate, 1e-7, 80, 300, None).unwrap();
    assert_eq!(report.stop, SqpStop::Converged, "{report:?}");
    for (actual, expected) in full.point().iter().zip(TRUTH) {
        assert!(
            (actual - expected).abs() < 2e-4 * (1.0 + expected.abs()),
            "recovered {actual}, expected {expected}"
        );
    }
    assert!(full.sample().f < initial * 1e-9);
    assert_eq!(
        split
            .try_run(&mut evaluate, 1e-7, 2, 300, None)
            .unwrap()
            .stop,
        SqpStop::IterationLimit
    );
    let mut fork = split.clone();
    let split_report = split.try_run(&mut evaluate, 1e-7, 1, 300, None).unwrap();
    let fork_report = fork.try_run(&mut evaluate, 1e-7, 1, 300, None).unwrap();
    assert_eq!(split_report.stop, fork_report.stop);
    same_accepted(&split, &fork);
    assert_eq!(split.evaluations(), fork.evaluations());
    assert_eq!(split.rejected_trials(), fork.rejected_trials());
    assert_eq!(
        split
            .try_run(&mut evaluate, 1e-7, 80, 300, None)
            .unwrap()
            .stop,
        SqpStop::Converged
    );
    same_accepted(&full, &split);
    assert_eq!(full.evaluations(), split.evaluations());
    assert_eq!(full.rejected_trials(), split.rejected_trials());
}

#[test]
fn cancellation_and_bad_observations_preserve_accepted_sample_and_charge_retries() {
    let mut data = PlateCalibration::synthetic().unwrap();
    let mut state = SqpState::try_new(
        &START,
        9,
        &mut |point| data.evaluate(point, &mut || false),
        None,
    )
    .unwrap();
    let before = state.clone();
    let mut polls = 0;
    let error = state
        .try_run(
            &mut |point| {
                data.evaluate(point, &mut || {
                    polls += 1;
                    polls >= 32
                })
            },
            1e-7,
            1,
            300,
            None,
        )
        .unwrap_err();
    assert!(
        matches!(error, SqpError::Evaluation(message) if message.to_ascii_lowercase().contains("cancel"))
    );
    assert!(polls >= 32);
    same_accepted(&before, &state);
    assert!(state.evaluations() > before.evaluations());

    let spent = state.evaluations();
    let target = data.targets[0];
    data.targets[0] = f64::NAN;
    assert!(matches!(
        state.try_run(
            &mut |point| data.evaluate(point, &mut || false),
            1e-7,
            1,
            300,
            None
        ),
        Err(SqpError::Evaluation(_))
    ));
    same_accepted(&before, &state);
    assert!(state.evaluations() > spent);
    data.targets[0] = target;

    let report = state
        .try_run(
            &mut |point| data.evaluate(point, &mut || false),
            1e-7,
            1,
            300,
            None,
        )
        .unwrap();
    assert!(matches!(
        report.stop,
        SqpStop::IterationLimit | SqpStop::Converged
    ));
    assert!(state.sample().f < before.sample().f);
    assert_eq!(state.iterations(), before.iterations() + 1);
    assert_eq!(
        state.sample(),
        &data
            .evaluate(state.point(), &mut || false)
            .unwrap()
            .unwrap()
    );
}

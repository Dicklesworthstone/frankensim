use super::*;
use crate::spectral::doppler::{VACUUM_LIGHT_SPEED_M_S as LIGHT, retarded::VacuumTrajectoryKnot};

fn dd(x: f64) -> Dd {
    Dd::from_f64(x)
}

fn uniform(x: f64, beta: f64, history: [f64; 2]) -> UniformVacuumWorldline {
    UniformVacuumWorldline::try_new(
        Dd::ZERO,
        Dd::ZERO,
        [x * LIGHT, 0.0, 0.0],
        [beta * LIGHT, 0.0, 0.0],
        history.map(dd),
    )
    .unwrap()
}

fn settings() -> RetardedExposureSettings {
    RetardedExposureSettings {
        max_ray_queries: 128,
        max_residual_m: 1e-6,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}

#[test]
fn g1_uniform_exposure_integrates_on_the_receiver_clock_without_a_second_doppler_factor() {
    let source = uniform(0.0, 0.6, [-10.0, 10.0]);
    let receiver = uniform(1.0, -0.2, [0.0, 1.0]);
    let mut calls = 0;
    let result = integrate_uniform(
        &source,
        &receiver,
        &[dd(0.0), dd(0.5)],
        settings(),
        |ray, _| {
            calls += 1;
            Ok::<_, ()>([
                1.0,
                ray.proper_time_rate(),
                (ray.emission_proper_time_s().to_f64() + 5.0).powi(2),
            ])
        },
        || false,
    )
    .unwrap();
    // Independent collinear law: D=sqrt(6), tau_emission=-2+D*tau_receiver.
    let d = 6.0_f64.sqrt();
    close(result.integrated[0], 0.5, 3e-15);
    close(result.integrated[1], 0.5 * d, 3e-15);
    close(result.integrated[2], 4.5 + 0.75 * d + d * d / 24.0, 5e-15);
    assert_eq!(result.receiver_proper_time_s, [dd(0.0), dd(0.5)]);
    assert_eq!(result.ray_queries, 6);
    assert_eq!(result.observable_evaluations, 4);
    assert_eq!(calls, 4);

    let common = integrate_uniform(
        &uniform(0.0, 0.6, [-10.0, 10.0]),
        &uniform(1.0, 0.6, [0.0, 1.0]),
        &[dd(0.0), dd(0.5)],
        settings(),
        |ray, _| Ok::<_, ()>([ray.proper_time_rate()]),
        || false,
    )
    .unwrap();
    close(common.integrated[0], 0.5, 3e-15);
}

#[test]
fn g1_piecewise_exposure_retains_each_detector_clock_rate_across_the_turn() {
    let knots = |values: &[(f64, f64)]| -> Vec<VacuumTrajectoryKnot> {
        values
            .iter()
            .map(|&(t, x)| VacuumTrajectoryKnot {
                coordinate_time_s: dd(t),
                position_m: [x * LIGHT, 0.0, 0.0],
            })
            .collect()
    };
    let source =
        PiecewiseVacuumWorldline::try_new(&knots(&[(-10.0, 0.0), (10.0, 0.0)]), dd(-10.0), 1)
            .unwrap();
    let receiver = PiecewiseVacuumWorldline::try_new(
        &knots(&[(0.0, 3.0), (2.0, 1.8), (4.0, 1.8)]),
        Dd::ZERO,
        2,
    )
    .unwrap();
    let end = receiver.proper_history_s()[1];
    let result = integrate_piecewise(
        &source,
        &receiver,
        &[dd(0.0), end - dd(2.0), end],
        settings(),
        |ray, _| Ok::<_, ()>([1.0, ray.proper_time_rate()]),
        || false,
    )
    .unwrap();
    // Receiver proper duration: 2*.8 + 2. Source emits from -3 to 2.2 s.
    close(result.integrated[0], 3.6, 4e-15);
    close(result.integrated[1], 5.2, 4e-15);
    assert_eq!(result.ray_queries, 10);
}

#[test]
fn g1_temporal_refinement_converges_without_claiming_an_error_certificate() {
    let source = uniform(0.0, 0.0, [-10.0, 10.0]);
    let receiver = uniform(1.0, 0.0, [0.0, 10.0]);
    let evaluate = |edges: &[Dd]| {
        integrate_uniform(
            &source,
            &receiver,
            edges,
            settings(),
            |ray, _| Ok::<_, ()>([ray.emission_proper_time_s().to_f64().exp()]),
            || false,
        )
        .unwrap()
    };
    let coarse = evaluate(&[dd(1.0), dd(2.0)]);
    let fine = evaluate(&[dd(1.0), dd(1.25), dd(1.5), dd(1.75), dd(2.0)]);
    let exact = 1.0_f64.exp_m1();
    let coarse_error = (coarse.integrated[0] - exact).abs();
    let fine_error = (fine.integrated[0] - exact).abs();
    assert!(coarse_error > 1e-11);
    assert!(fine_error < 2e-14 && fine_error < 1e-4 * coarse_error);
    assert_eq!(
        fine,
        evaluate(&[dd(1.0), dd(1.25), dd(1.5), dd(1.75), dd(2.0)])
    );
}

#[test]
fn g3_full_history_and_work_admission_precede_observable_callbacks() {
    let source = uniform(0.0, 0.0, [0.0, 1.0]);
    let receiver = uniform(1.0, 0.0, [0.0, 10.0]);
    let mut calls = 0;
    // All four interior nodes would be within the source history; the final
    // unsampled tail is not. Endpoint preflight must refuse before reading it.
    let result = integrate_uniform(
        &source,
        &receiver,
        &[dd(1.0), dd(2.01)],
        settings(),
        |_, _| {
            calls += 1;
            Ok::<_, ()>([1.0])
        },
        || false,
    );
    assert_eq!(
        result,
        Err(RetardedExposureError::Geometry(
            RetardedVacuumError::OutsideHistory
        ))
    );
    assert_eq!(calls, 0);
    let mut budget = settings();
    budget.max_ray_queries = 5;
    let result = integrate_uniform(
        &source,
        &receiver,
        &[dd(1.0), dd(2.0)],
        budget,
        |_, _| {
            calls += 1;
            Ok::<_, ()>([1.0])
        },
        || false,
    );
    assert_eq!(result, Err(RetardedExposureError::RayBudget));
    assert_eq!(calls, 0);
    budget.max_ray_queries = 6;
    assert!(
        integrate_uniform(
            &source,
            &receiver,
            &[dd(1.0), dd(2.0)],
            budget,
            |_, _| Ok::<_, ()>([1.0]),
            || false
        )
        .is_ok()
    );
    for edges in [
        &[][..],
        &[dd(1.0)][..],
        &[dd(1.0), dd(1.0)][..],
        &[dd(2.0), dd(1.0)][..],
    ] {
        assert_eq!(
            integrate_uniform(
                &source,
                &receiver,
                edges,
                settings(),
                |_, _| Ok::<_, ()>([1.0]),
                || false
            ),
            Err(RetardedExposureError::InvalidPanels)
        );
    }
}

#[test]
fn g4_cancellation_and_invalid_rates_never_publish_partial_exposure() {
    let source = uniform(0.0, 0.0, [-10.0, 10.0]);
    let receiver = uniform(1.0, 0.0, [0.0, 10.0]);
    let mut preflight_checks = 0;
    let preflight = integrate_uniform(
        &source,
        &receiver,
        &[dd(0.0), dd(1.0), dd(2.0), dd(2.0)],
        settings(),
        |_, _| Ok::<_, ()>([1.0]),
        || {
            preflight_checks += 1;
            preflight_checks >= 2
        },
    );
    // A cancelled preflight must stop before reading the later invalid panel.
    assert_eq!(preflight, Err(RetardedExposureError::Cancelled));
    assert_eq!(preflight_checks, 2);
    let mut checkpoints = 0;
    let mut samples = 0;
    let result = integrate_uniform(
        &source,
        &receiver,
        &[dd(1.0), dd(2.0)],
        settings(),
        |_, _| {
            samples += 1;
            Ok::<_, ()>([1.0])
        },
        || {
            checkpoints += 1;
            checkpoints == 6
        },
    );
    assert_eq!(result, Err(RetardedExposureError::Cancelled));
    assert_eq!(samples, 2);
    checkpoints = 0;
    samples = 0;
    let in_callback = integrate_uniform(
        &source,
        &receiver,
        &[dd(1.0), dd(2.0)],
        settings(),
        |_, checkpoint| {
            samples += 1;
            if checkpoint() {
                Err("cancelled in observable")
            } else {
                Ok([1.0])
            }
        },
        || {
            checkpoints += 1;
            checkpoints >= 5
        },
    );
    assert_eq!(
        in_callback,
        Err(RetardedExposureError::Observable("cancelled in observable"))
    );
    assert_eq!(samples, 1);
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            integrate_uniform(
                &source,
                &receiver,
                &[dd(1.0), dd(2.0)],
                settings(),
                |_, _| Ok::<_, ()>([invalid]),
                || false
            ),
            Err(RetardedExposureError::InvalidRate)
        );
    }
    assert_eq!(
        integrate_uniform(
            &source,
            &receiver,
            &[dd(1.0), dd(2.0)],
            settings(),
            |_, _| Err::<[f64; 1], _>("source unavailable"),
            || false
        ),
        Err(RetardedExposureError::Observable("source unavailable"))
    );
}

#[test]
fn g3_large_receiver_epochs_keep_the_actual_small_exposure_duration() {
    let epoch = dd(1e20);
    let source = uniform(0.0, 0.0, [-10.0, 10.0]);
    let receiver = UniformVacuumWorldline::try_new(
        Dd::ZERO,
        epoch,
        [LIGHT, 0.0, 0.0],
        [0.0; 3],
        [epoch, epoch + dd(1.0)],
    )
    .unwrap();
    let result = integrate_uniform(
        &source,
        &receiver,
        &[epoch, epoch + dd(0.001)],
        settings(),
        |_, _| Ok::<_, ()>([1.0]),
        || false,
    )
    .unwrap();
    close(result.integrated[0], 0.001, 3e-15);
}

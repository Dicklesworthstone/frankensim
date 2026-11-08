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

fn piecewise(values: &[(f64, f64)], proper_epoch: Dd) -> PiecewiseVacuumWorldline {
    let knots: Vec<_> = values
        .iter()
        .map(|&(time, x)| VacuumTrajectoryKnot {
            coordinate_time_s: dd(time),
            position_m: [x * LIGHT, 0.0, 0.0],
        })
        .collect();
    PiecewiseVacuumWorldline::try_new(&knots, proper_epoch, values.len() - 1).unwrap()
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
fn g1_motion_alignment_splits_at_the_received_source_turn_not_its_emission_time() {
    let source = piecewise(&[(0.0, 0.0), (1.0, 0.0), (4.0, 1.8)], Dd::ZERO);
    let receiver = piecewise(&[(3.0, 3.0), (5.0, 3.0)], dd(3.0));
    let edges = [dd(3.8), dd(4.8)];
    let rates = |ray: &RetardedVacuumRay, _: &mut dyn FnMut() -> bool| {
        Ok::<_, ()>([1.0, ray.proper_time_rate()])
    };
    let original =
        integrate_piecewise(&source, &receiver, &edges, settings(), rates, || false).unwrap();
    let aligned =
        integrate_piecewise_with_knots(&source, &receiver, &edges, settings(), rates, || false)
            .unwrap();
    // Turn at emission t=1 arrives at t=4. D=1 for .2 s, then D=2 for .8 s.
    close(aligned.integrated[0], 1.0, 3e-15);
    close(aligned.integrated[1], 1.8, 5e-15);
    assert!((original.integrated[1] - 1.8).abs() > 0.02);
    assert_eq!(original.ray_queries, 6);
    assert_eq!(aligned.ray_queries, 13); // two planning ends + one arrival + 2 + 8
    assert_eq!(aligned.observable_evaluations, 8);
}

#[test]
fn g1_motion_alignment_uses_accumulated_detector_proper_time_at_a_large_epoch() {
    let source = piecewise(&[(-10.0, 0.0), (10.0, 0.0)], dd(-10.0));
    let epoch = dd(1e20);
    let receiver = piecewise(&[(0.0, 3.0), (2.0, 1.8), (4.0, 1.8)], epoch);
    let result = integrate_piecewise_with_knots(
        &source,
        &receiver,
        &receiver.proper_history_s(),
        settings(),
        |ray, _| Ok::<_, ()>([1.0, ray.proper_time_rate()]),
        || false,
    )
    .unwrap();
    // The detector's turn is at tau=epoch+1.6, not coordinate t=2.
    close(result.integrated[0], 3.6, 4e-15);
    close(result.integrated[1], 5.2, 4e-15);
    assert_eq!(result.ray_queries, 12);
    assert_eq!(result.observable_evaluations, 8);
}

#[test]
fn g3_motion_alignment_filters_irrelevant_knots_and_deduplicates_one_physical_edge() {
    let source = piecewise(
        &[
            (-10.0, 0.0),
            (-5.0, 0.0),
            (0.0, 0.0),
            (1.0, 0.0),
            (2.0, 0.0),
            (5.0, 0.0),
            (10.0, 0.0),
        ],
        dd(-10.0),
    );
    let receiver = piecewise(&[(3.5, 3.0), (4.0, 3.0), (4.5, 3.0)], dd(3.5));
    let mut budget = settings();
    budget.max_ray_queries = 13;
    let result = integrate_piecewise_with_knots(
        &source,
        &receiver,
        &[dd(3.5), dd(4.0), dd(4.5)],
        budget,
        |_, _| Ok::<_, ()>([1.0]),
        || false,
    )
    .unwrap();
    // All other source-knot arrivals are outside the receiver history. The
    // authored edge, detector knot and relevant source arrival coincide at 4.
    close(result.integrated[0], 1.0, 3e-15);
    assert_eq!(result.ray_queries, 13);
    assert_eq!(result.observable_evaluations, 8);
    budget.max_ray_queries = 12;
    let mut calls = 0;
    let refused = integrate_piecewise_with_knots(
        &source,
        &receiver,
        &[dd(3.5), dd(4.5)],
        budget,
        |_, _| {
            calls += 1;
            Ok::<_, ()>([1.0])
        },
        || false,
    );
    assert_eq!(refused, Err(RetardedExposureError::RayBudget));
    assert_eq!(calls, 0);
    budget.max_ray_queries = 8;
    let clear = integrate_piecewise_with_knots(
        &source,
        &receiver,
        &[dd(3.6), dd(3.9)],
        budget,
        |_, _| Ok::<_, ()>([1.0]),
        || false,
    )
    .unwrap();
    close(clear.integrated[0], 0.3, 3e-15);
    assert_eq!(clear.ray_queries, 8); // no out-of-exposure arrival queries
}

#[test]
fn g4_motion_planning_cancels_and_refuses_missing_history_before_callbacks() {
    let source = piecewise(
        &[
            (-10.0, 0.0),
            (0.0, 0.0),
            (1.0, 0.0),
            (2.0, 0.0),
            (10.0, 0.0),
        ],
        dd(-10.0),
    );
    let receiver = piecewise(&[(3.5, 3.0), (4.0, 3.0), (4.5, 3.0)], dd(3.5));
    for stop in [4, 11] {
        let mut checks = 0;
        let mut calls = 0;
        let cancelled = integrate_piecewise_with_knots(
            &source,
            &receiver,
            &[dd(3.5), dd(4.5)],
            settings(),
            |_, _| {
                calls += 1;
                Ok::<_, ()>([1.0])
            },
            || {
                checks += 1;
                checks == stop
            },
        );
        assert_eq!(cancelled, Err(RetardedExposureError::Cancelled));
        assert_eq!(checks, stop);
        assert_eq!(calls, 0);
    }
    let mut calls = 0;
    let missing = integrate_piecewise_with_knots(
        &source,
        &receiver,
        &[dd(3.5), dd(4.51)],
        settings(),
        |_, _| {
            calls += 1;
            Ok::<_, ()>([1.0])
        },
        || false,
    );
    assert_eq!(
        missing,
        Err(RetardedExposureError::Geometry(
            RetardedVacuumError::OutsideHistory
        ))
    );
    assert_eq!(calls, 0);
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

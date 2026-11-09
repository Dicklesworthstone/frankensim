//! Source boundary clocks, signed affine prescriptions and explicit separation.
use super::{close, controls, directions, driver, source};
use fs_contact::normal_patch::adhesive::mixed_mode::combined::{
    JkrCombinedBoundary,
    loading::{JkrCombinedLoading, JkrCombinedLoadingControls, JkrCombinedLoadingError},
};
use fs_exec::CancelGate;

#[test]
fn both_loading_directions_retain_the_original_step_and_locate_the_first_physical_boundary() {
    let source = source();
    let gate = CancelGate::new_clock_free();
    let (a, direction) = directions(&source, 4.0, 8.0);
    let delta = -a * a / source.normal_source().effective_radius_m();
    let loading = driver(&source, delta, direction, &gate);
    let limit = loading.limit().unwrap().load_factor;
    let initial = loading.initial_state(0.0, &gate).unwrap();
    for sign in [-1.0, 1.0] {
        let requested = sign * 2.0 * limit;
        let expected = loading
            .advance_to_load_factor(initial, requested, &gate)
            .unwrap();
        let ramp = loading
            .advance_linear_ramp(initial, 2.0, 6.0, requested, 0.0, &gate)
            .unwrap();
        assert_eq!(ramp.step(), &expected);
        assert_eq!(ramp.previous_state(), initial);
        assert_eq!(ramp.state(), expected.state());
        assert_eq!(ramp.initial_time_s(), 2.0);
        assert_eq!(ramp.requested_final_time_s(), 6.0);
        assert_eq!(ramp.accepted_final_time_s(), 4.0);
        assert_eq!(ramp.boundary_time_s(), Some(4.0));
        assert_eq!(ramp.untraversed_time_s(), 2.0);
        assert_eq!(ramp.clock_residual_s(), 0.0);
        assert_eq!(ramp.maximum_clock_error_s(), 0.0);
        assert_eq!(ramp.step().report().released_energy_j, 0.0);
        assert_eq!(
            ramp.state().contact().unwrap().boundary,
            Some(JkrCombinedBoundary::RadialFold)
        );
        for (time, fraction) in [(2.0, 0.0), (3.0, 0.5), (4.0, 1.0), (5.0, 1.5), (6.0, 2.0)] {
            close(
                ramp.load_factor_at_time(time).unwrap(),
                sign * fraction * limit,
                3e-15,
            );
        }
        assert_eq!(
            ramp.load_factor_at_time(4.0).unwrap(),
            ramp.state().load_factor()
        );
        assert_eq!(ramp.load_factor_at_time(6.0).unwrap(), requested);
        assert!(
            loading
                .detach_at_fold(ramp.state(), &gate)
                .unwrap()
                .report()
                .released_energy_j
                > 0.0
        );
    }
    let exact_arrival = loading
        .advance_linear_ramp(initial, 2.0, 6.0, limit, 0.0, &gate)
        .unwrap();
    assert_eq!(exact_arrival.boundary_time_s(), Some(6.0));
    assert_eq!(exact_arrival.untraversed_time_s(), 0.0);
    assert_eq!(exact_arrival.load_factor_at_time(6.0).unwrap(), limit);
    assert_eq!(loading.initial_state(0.0, &gate).unwrap(), initial);
}

#[test]
fn reversal_crosses_zero_on_the_prescribed_clock_without_losing_the_radius_extremum() {
    let source = source();
    let gate = CancelGate::new_clock_free();
    let (a, direction) = directions(&source, 4.0, 8.0);
    let delta = -a * a / source.normal_source().effective_radius_m();
    let loading = driver(&source, delta, direction, &gate);
    let limit = loading.limit().unwrap().load_factor;
    let initial = loading.initial_state(0.5 * limit, &gate).unwrap();
    let requested = -2.0 * limit;
    let ramp = loading
        .advance_linear_ramp(initial, 7.0, 12.0, requested, 1e-13, &gate)
        .unwrap();
    assert_eq!(ramp.state().load_factor(), -limit);
    close(ramp.accepted_final_time_s(), 10.0, 4e-15);
    assert_eq!(ramp.boundary_time_s(), Some(ramp.accepted_final_time_s()));
    assert!(ramp.load_factor_at_time(8.0).unwrap().abs() <= 8.0 * f64::EPSILON * limit);
    assert_eq!(
        ramp.step().report().maximum_contact_radius_m,
        loading
            .initial_state(0.0, &gate)
            .unwrap()
            .contact()
            .unwrap()
            .contact_radius_m
            .next_up()
    );
    assert_eq!(ramp.step().report().evaluations, 1);
    assert_eq!(ramp.step().report().released_energy_j, 0.0);
    // The timestamp follows the actual retained endpoints, including their SI
    // floating-point representation, rather than a separately rounded rate.
    let fraction =
        (ramp.state().load_factor() - initial.load_factor()) / (requested - initial.load_factor());
    close(ramp.accepted_final_time_s(), 7.0 + 5.0 * fraction, 4e-15);
}

#[test]
fn held_and_outward_boundary_requests_and_detached_tails_keep_distinct_physical_meanings() {
    let source = source();
    let gate = CancelGate::new_clock_free();
    let (a, direction) = directions(&source, 4.0, 8.0);
    let delta = -a * a / source.normal_source().effective_radius_m();
    let loading = driver(&source, delta, direction, &gate);
    let limit = loading.limit().unwrap().load_factor;
    let fold = loading.initial_state(limit, &gate).unwrap();
    let outward = loading
        .advance_linear_ramp(fold, 100.0, 102.0, 2.0 * limit, 0.0, &gate)
        .unwrap();
    assert_eq!(outward.accepted_final_time_s(), 100.0);
    assert_eq!(outward.boundary_time_s(), Some(100.0));
    assert_eq!(outward.untraversed_time_s(), 2.0);
    assert_eq!(outward.state().contact(), fold.contact());
    assert_eq!(outward.step().report().evaluations, 0);
    assert_eq!(outward.step().report().released_energy_j, 0.0);
    let held = loading
        .advance_linear_ramp(fold, 100.0, 102.0, limit, 0.0, &gate)
        .unwrap();
    assert_eq!(held.accepted_final_time_s(), 102.0);
    assert_eq!(held.boundary_time_s(), None);
    assert_eq!(held.untraversed_time_s(), 0.0);
    assert_eq!(held.state().contact(), fold.contact());
    assert_eq!(held.step().report().prescribed_work_j, 0.0);
    assert_eq!(held.load_factor_at_time(101.0).unwrap(), limit);
    let detached = loading
        .detach_at_fold(outward.state(), &gate)
        .unwrap()
        .state();
    let free = loading
        .advance_linear_ramp(detached, 100.0, 102.0, -2.0 * limit, 0.0, &gate)
        .unwrap();
    assert_eq!(free.accepted_final_time_s(), 102.0);
    assert_eq!(free.boundary_time_s(), None);
    assert_eq!(free.state().contact(), None);
    assert_eq!(
        free.state().released_energy_j(),
        detached.released_energy_j()
    );
    assert_eq!(free.step().report().evaluations, 0);
    assert_eq!(free.step().report().prescribed_work_j, 0.0);
    assert_eq!(free.step().report().released_energy_j, 0.0);
    assert_eq!(free.state().steps(), detached.steps() + 1);

    let (a, direction) = directions(&source, 8.0, 2.0);
    let positive = a * a / source.normal_source().effective_radius_m();
    let opening = driver(&source, positive, direction, &gate);
    let limit = opening.limit().unwrap().load_factor;
    let initial = opening.initial_state(0.0, &gate).unwrap();
    let ramp = opening
        .advance_linear_ramp(initial, 10.0, 12.0, 2.0 * limit, 0.0, &gate)
        .unwrap();
    assert_eq!(ramp.boundary_time_s(), Some(11.0));
    assert_eq!(
        ramp.state().contact().unwrap().boundary,
        Some(JkrCombinedBoundary::OpeningEdge)
    );
    assert!(matches!(
        opening.detach_at_fold(ramp.state(), &gate),
        Err(JkrCombinedLoadingError::NotAtSeparatingFold)
    ));
}

#[test]
fn collapsed_boundary_clocks_and_excess_roundoff_refuse_without_moving_or_releasing_history() {
    let source = source();
    let gate = CancelGate::new_clock_free();
    let (a, direction) = directions(&source, 4.0, 8.0);
    let delta = -a * a / source.normal_source().effective_radius_m();
    let loading = driver(&source, delta, direction, &gate);
    let limit = loading.limit().unwrap().load_factor;
    let initial = loading.initial_state(0.0, &gate).unwrap();
    let epoch = 2.0_f64.powi(53);
    assert!(matches!(
        loading.advance_linear_ramp(initial, epoch, epoch.next_up(), 2.0 * limit, 10.0, &gate),
        Err(JkrCombinedLoadingError::Invalid(_))
    ));
    assert!(matches!(
        loading.advance_linear_ramp(initial, 0.0, 0.0_f64.next_up(), 2.0 * limit, 1.0, &gate),
        Err(JkrCombinedLoadingError::Invalid(_))
    ));
    let requested = 3.0 * limit;
    assert!(matches!(
        loading.advance_linear_ramp(initial, 1.0, 2.0, requested, 0.0, &gate),
        Err(JkrCombinedLoadingError::ClockRoundoff { residual_s, tolerance_s })
            if residual_s != 0.0 && residual_s.is_finite() && tolerance_s == 0.0
    ));
    let tolerance = 8.0 * f64::EPSILON;
    let admitted = loading
        .advance_linear_ramp(initial, 1.0, 2.0, requested, tolerance, &gate)
        .unwrap();
    assert!(admitted.clock_residual_s().abs() <= tolerance);
    assert_ne!(admitted.clock_residual_s(), 0.0);
    assert!(admitted.accepted_final_time_s() > 1.0 && admitted.accepted_final_time_s() < 2.0);
    assert_eq!(
        admitted
            .load_factor_at_time(admitted.accepted_final_time_s())
            .unwrap(),
        limit
    );
    for (start, end, tolerance) in [
        (-1.0, 2.0, 0.0),
        (1.0, 1.0, 0.0),
        (2.0, 1.0, 0.0),
        (f64::NAN, 2.0, 0.0),
        (0.0, f64::INFINITY, 0.0),
        (0.0, 2.0, -1.0),
        (0.0, 2.0, f64::NAN),
        (0.0, 2.0, f64::INFINITY),
    ] {
        assert!(matches!(
            loading.advance_linear_ramp(initial, start, end, requested, tolerance, &gate),
            Err(JkrCombinedLoadingError::Invalid(_))
        ));
    }
    for bad in [f64::NAN, f64::INFINITY, 0.0, 3.0] {
        assert!(admitted.load_factor_at_time(bad).is_err());
    }
    assert_eq!(loading.initial_state(0.0, &gate).unwrap(), initial);
}

#[test]
fn unbounded_twist_family_uses_one_original_query_and_rejects_foreign_or_cancelled_history() {
    let source = source();
    let gate = CancelGate::new_clock_free();
    let (_, direction) = directions(&source, 4.0, 0.0);
    let loading = JkrCombinedLoading::new(
        &source,
        0.0,
        [0.0; 2],
        direction[2],
        JkrCombinedLoadingControls {
            maximum_evaluations: 1,
            ..controls()
        },
        &gate,
    )
    .unwrap();
    assert_eq!(loading.limit(), None);
    let initial = loading.initial_state(0.2, &gate).unwrap();
    let ramp = loading
        .advance_linear_ramp(initial, 0.0, 0.25, 0.3, 0.0, &gate)
        .unwrap();
    assert_eq!(ramp.step().report().evaluations, 1);
    assert_eq!(ramp.accepted_final_time_s(), 0.25);
    assert_eq!(ramp.boundary_time_s(), None);
    assert_eq!(ramp.untraversed_time_s(), 0.0);
    close(ramp.load_factor_at_time(0.125).unwrap(), 0.25, 3e-15);
    let foreign = driver(&source, 0.0, [0.0, 0.0, -direction[2]], &gate);
    assert!(matches!(
        foreign.advance_linear_ramp(initial, 0.0, 1.0, 0.3, 0.0, &gate),
        Err(JkrCombinedLoadingError::StateMismatch)
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        loading.advance_linear_ramp(initial, 0.0, 1.0, 0.3, 0.0, &cancelled),
        Err(JkrCombinedLoadingError::Cancelled)
    ));
    assert_eq!(loading.initial_state(0.2, &gate).unwrap(), initial);
}

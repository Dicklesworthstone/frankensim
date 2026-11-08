//! First twist boundary, reversible work and explicitly retained snap release.
use super::*;
use fs_contact::normal_patch::adhesive::mixed_mode::torsion::loading::*;

fn controls() -> JkrTorsionalLoadingControls {
    JkrTorsionalLoadingControls {
        maximum_evaluations: 2,
        energy_absolute_tolerance_j: 0.0,
        energy_relative_tolerance: 2e-12,
    }
}

#[test]
fn prescribed_twist_reaches_both_radial_folds_and_releases_original_energy_once() {
    let source = source();
    let gate = CancelGate::new_clock_free();
    let delta = 0.6 * source.normal_source().minimum_attached_indentation_m();
    let loading = JkrTorsionalLoading::new(&source, delta, controls(), &gate).unwrap();
    let limit = loading.limit().unwrap();
    assert_eq!(limit.kind, JkrTorsionalLimitKind::RadialFold);
    let mut releases = [0.0; 2];
    for (index, sign) in [-1.0, 1.0].into_iter().enumerate() {
        let initial = loading.initial_state(0.0, &gate).unwrap();
        let requested = sign * 2.0 * limit.twist_angle_rad;
        let step = loading.advance_to_twist(initial, requested, &gate).unwrap();
        let report = step.report();
        let at_fold = step.state();
        let fold = at_fold.contact().unwrap();
        assert_eq!(step.previous_state(), initial);
        assert_eq!(at_fold.twist_angle_rad(), sign * limit.twist_angle_rad);
        assert_eq!(report.untraversed_twist_rad, sign * limit.twist_angle_rad);
        assert_eq!(report.boundary, Some(JkrTorsionalLimitKind::RadialFold));
        assert_eq!(fold.displacement_radius_curvature_n_per_m, 0.0);
        assert!(fold.coupled_stiffness.is_none());
        assert_eq!(report.released_energy_j, 0.0);
        assert!(report.prescribed_work_j > 0.0);
        assert!(report.energy_residual_j.abs() <= report.energy_tolerance_j);
        // Independent fold elimination: no heat fraction or normal-only snap.
        let normal = source.normal_source();
        let e = normal.normal_state().elastic().reduced_modulus_pa();
        let a = fold.contact_radius_m;
        let hertz = a * a / normal.effective_radius_m();
        let expected_release = e * a * (delta * delta / 3.0 + hertz * hertz / 5.0);
        let detached = loading.detach_at_fold(at_fold, &gate).unwrap();
        releases[index] = detached.report().released_energy_j;
        close(releases[index], expected_release, 2e-12);
        assert_eq!(detached.report().prescribed_work_j, 0.0);
        assert_eq!(detached.report().energy_residual_j, 0.0);
        assert_eq!(detached.state().contact(), None);
        assert_eq!(
            detached.state().twist_angle_rad(),
            at_fold.twist_angle_rad()
        );
        assert_eq!(detached.state().released_energy_j(), releases[index]);
        assert_eq!(detached.state().stored_energy_j(), 0.0);
        close(
            initial.stored_energy_j() + report.prescribed_work_j,
            releases[index],
            2e-12,
        );
        assert!(matches!(
            loading.detach_at_fold(detached.state(), &gate),
            Err(JkrTorsionalLoadingError::NotAtRadialFold)
        ));
        let free = loading
            .advance_to_twist(detached.state(), -requested, &gate)
            .unwrap();
        assert_eq!(free.state().contact(), None);
        assert_eq!(free.report().prescribed_work_j, 0.0);
        assert_eq!(free.report().released_energy_j, 0.0);
        assert_eq!(free.report().untraversed_twist_rad, 0.0);
        assert_eq!(free.state().released_energy_j(), releases[index]);
        assert_eq!(free.state().steps(), initial.steps() + 3);
    }
    assert_eq!(releases[0], releases[1]);
}

#[test]
fn reversible_twist_work_matches_independent_torque_quadrature_and_reversal_radius() {
    let source = source();
    let gate = CancelGate::new_clock_free();
    let delta = 0.5 * source.normal_source().minimum_attached_indentation_m();
    let loading = JkrTorsionalLoading::new(&source, delta, controls(), &gate).unwrap();
    let phi = loading.limit().unwrap().twist_angle_rad;
    let from = loading.initial_state(0.15 * phi, &gate).unwrap();
    let forward = loading.advance_to_twist(from, 0.7 * phi, &gate).unwrap();
    let end = forward.state();
    let increment = end.twist_angle_rad() - from.twist_angle_rad();
    let intervals = 512;
    let mut torque_work = 0.0;
    for i in 0..=intervals {
        let angle = from.twist_angle_rad() + increment * (i as f64 / intervals as f64);
        let point = source
            .evaluate_torsion_displacement_controlled(delta, angle)
            .unwrap();
        let weight = if i == 0 || i == intervals {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        torque_work += weight * point.torsional_moment_n_m;
        assert!(point.contact_radius_m >= forward.report().minimum_contact_radius_m);
        assert!(point.contact_radius_m <= forward.report().maximum_contact_radius_m);
    }
    torque_work *= increment / (3.0 * intervals as f64);
    close(forward.report().prescribed_work_j, torque_work, 2e-9);
    let reverse = loading
        .advance_to_twist(end, from.twist_angle_rad(), &gate)
        .unwrap();
    assert_eq!(
        reverse.report().prescribed_work_j,
        -forward.report().prescribed_work_j
    );
    assert_eq!(reverse.state().contact(), from.contact());
    assert_eq!(reverse.state().released_energy_j(), 0.0);
    let crossing = loading.advance_to_twist(end, -0.4 * phi, &gate).unwrap();
    let zero = source
        .evaluate_torsion_displacement_controlled(delta, 0.0)
        .unwrap();
    assert_eq!(crossing.report().evaluations, 2);
    assert_eq!(
        crossing.report().maximum_contact_radius_m,
        zero.contact_radius_m.next_up()
    );
    assert!(zero.contact_radius_m > end.contact().unwrap().contact_radius_m);
    close(
        crossing.report().prescribed_work_j,
        crossing.state().stored_energy_j() - end.stored_energy_j(),
        2e-12,
    );
}

#[test]
fn opening_edge_is_an_applicability_boundary_and_zero_approach_has_no_invented_fold() {
    let source = source();
    let gate = CancelGate::new_clock_free();
    let (_, delta, phi) = oracle(&source);
    assert!(delta > 0.0);
    let loading = JkrTorsionalLoading::new(&source, delta, controls(), &gate).unwrap();
    let limit = loading.limit().unwrap();
    assert_eq!(limit.kind, JkrTorsionalLimitKind::OpeningEdge);
    let initial = loading.initial_state(0.0, &gate).unwrap();
    let stopped = loading
        .advance_to_twist(initial, 2.0 * limit.twist_angle_rad, &gate)
        .unwrap();
    assert_eq!(
        stopped.report().boundary,
        Some(JkrTorsionalLimitKind::OpeningEdge)
    );
    assert_eq!(stopped.report().released_energy_j, 0.0);
    assert!(matches!(
        loading.detach_at_fold(stopped.state(), &gate),
        Err(JkrTorsionalLoadingError::NotAtRadialFold)
    ));
    let unloaded = loading
        .advance_to_twist(stopped.state(), 0.0, &gate)
        .unwrap();
    assert_eq!(unloaded.state().contact(), initial.contact());
    assert_eq!(
        unloaded.report().prescribed_work_j,
        -stopped.report().prescribed_work_j
    );
    let zero = JkrTorsionalLoading::new(&source, 0.0, controls(), &gate).unwrap();
    assert_eq!(zero.limit(), None);
    let start = zero.initial_state(0.0, &gate).unwrap();
    let moved = zero.advance_to_twist(start, phi, &gate).unwrap();
    assert_eq!(moved.report().boundary, None);
    assert_eq!(moved.report().untraversed_twist_rad, 0.0);
    assert!(
        moved.state().contact().unwrap().contact_area_m2 < start.contact().unwrap().contact_area_m2
    );
}

#[test]
fn tiny_twist_work_is_retained_and_fold_work_needs_no_singular_tangent() {
    let source = source();
    let gate = CancelGate::new_clock_free();
    let delta = 0.7 * source.normal_source().minimum_attached_indentation_m();
    let loading = JkrTorsionalLoading::new(&source, delta, controls(), &gate).unwrap();
    let limit = loading.limit().unwrap().twist_angle_rad;
    let a = 0.3 * limit;
    let initial = loading.initial_state(a, &gate).unwrap();
    let tiny = loading
        .advance_to_twist(initial, a.next_up(), &gate)
        .unwrap();
    assert!(tiny.report().prescribed_work_j > 0.0);
    close(
        tiny.report().prescribed_work_j / (a.next_up() - a),
        initial.contact().unwrap().torsional_moment_n_m,
        3e-12,
    );
    let near_fold = loading.initial_state(limit.next_down(), &gate).unwrap();
    let final_step = loading.advance_to_twist(near_fold, limit, &gate).unwrap();
    assert!(final_step.report().prescribed_work_j > 0.0);
    assert_eq!(
        final_step.report().boundary,
        Some(JkrTorsionalLimitKind::RadialFold)
    );
    assert!(final_step.report().energy_residual_j.abs() <= final_step.report().energy_tolerance_j);
    let returned = loading
        .advance_to_twist(final_step.state(), limit.next_down(), &gate)
        .unwrap();
    assert_eq!(
        returned.report().prescribed_work_j,
        -final_step.report().prescribed_work_j
    );
}

#[test]
fn bounded_queries_foreign_sources_and_cancellation_leave_exact_retry_available() {
    let source = source();
    let gate = CancelGate::new_clock_free();
    let delta = 0.5 * source.normal_source().minimum_attached_indentation_m();
    let loading = JkrTorsionalLoading::new(&source, delta, controls(), &gate).unwrap();
    let phi = loading.limit().unwrap().twist_angle_rad;
    let previous = loading.initial_state(0.3 * phi, &gate).unwrap();
    let expected = loading
        .advance_to_twist(previous, -0.3 * phi, &gate)
        .unwrap();
    let tight = JkrTorsionalLoading::new(
        &source,
        delta,
        JkrTorsionalLoadingControls {
            maximum_evaluations: 1,
            ..controls()
        },
        &gate,
    )
    .unwrap();
    assert!(matches!(
        tight.advance_to_twist(previous, -0.3 * phi, &gate),
        Err(JkrTorsionalLoadingError::Budget)
    ));
    let other_approach = JkrTorsionalLoading::new(&source, delta * 0.9, controls(), &gate).unwrap();
    assert!(matches!(
        other_approach.advance_to_twist(previous, 0.4 * phi, &gate),
        Err(JkrTorsionalLoadingError::StateMismatch)
    ));
    let other_source = SourceJkrMixedMode::new(
        source.normal_source().clone(),
        JkrMixedModeEnvelope {
            maximum_tangential_displacement_to_radius: 0.09,
            ..source.envelope()
        },
    )
    .unwrap();
    let foreign = JkrTorsionalLoading::new(&other_source, delta, controls(), &gate).unwrap();
    assert!(matches!(
        foreign.advance_to_twist(previous, 0.4 * phi, &gate),
        Err(JkrTorsionalLoadingError::StateMismatch)
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        loading.advance_to_twist(previous, phi, &cancelled),
        Err(JkrTorsionalLoadingError::Cancelled)
    ));
    assert!(matches!(
        loading.detach_at_fold(previous, &gate),
        Err(JkrTorsionalLoadingError::NotAtRadialFold)
    ));
    for bad in [f64::NAN, f64::INFINITY] {
        assert!(loading.advance_to_twist(previous, bad, &gate).is_err());
    }
    assert_eq!(
        loading
            .advance_to_twist(previous, -0.3 * phi, &gate)
            .unwrap(),
        expected
    );
    assert_eq!(loading.initial_state(0.3 * phi, &gate).unwrap(), previous);
}

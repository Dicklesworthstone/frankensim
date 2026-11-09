//! Combined first-boundary physics, reversible work and explicit release.
#[path = "loading/temperature.rs"]
mod temperature;
#[path = "loading/timed.rs"]
mod timed;

use super::{close, source};
use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::mixed_mode::{
    JkrMixedModeEnvelope, SourceJkrMixedMode,
    combined::{JkrCombinedBoundary, JkrCombinedResponse, loading::*},
    torsion::JkrTorsionalLimitKind,
};
use fs_exec::CancelGate;

fn controls() -> JkrCombinedLoadingControls {
    JkrCombinedLoadingControls {
        maximum_iterations: 128,
        maximum_evaluations: 192,
        energy_absolute_tolerance_j: 0.0,
        energy_relative_tolerance: 3e-12,
    }
}

// Independent dimensional inputs producing the specified b,t in
// F=(y²-d)²+lambda²(b*y²+t)-16*y. No combined root/limit API enters.
fn directions(source: &SourceJkrMixedMode, b: f64, t: f64) -> (f64, [f64; 3]) {
    let normal = source.normal_source();
    let e = normal.normal_state().elastic().reduced_modulus_pa();
    let w = normal.normal_state().adhesion_energy_j_per_m2();
    let r = normal.effective_radius_m();
    let a = (PI * w * r * r / (8.0 * e)).cbrt();
    let depth = a * a / r;
    let stationary_scale = e * depth * depth;
    let shear = (t * stationary_scale / (4.0 * source.reduced_shear_modulus_pa())).sqrt();
    let twist =
        (b * stationary_scale / (8.0 * source.torsional_reduced_shear_modulus_pa() * a * a)).sqrt();
    (a, [0.6 * shear, 0.8 * shear, twist])
}

fn driver<'a>(
    source: &'a SourceJkrMixedMode,
    delta: f64,
    direction: [f64; 3],
    gate: &CancelGate,
) -> JkrCombinedLoading<'a> {
    JkrCombinedLoading::new(
        source,
        delta,
        [direction[0], direction[1]],
        direction[2],
        controls(),
        gate,
    )
    .unwrap()
}

fn potential(source: &SourceJkrMixedMode, point: JkrCombinedResponse) -> [f64; 4] {
    let normal = source.normal_source();
    let e = normal.normal_state().elastic().reduced_modulus_pa();
    let w = normal.normal_state().adhesion_energy_j_per_m2();
    let r = normal.effective_radius_m();
    let a = point.contact_radius_m;
    let q = point.coordinates;
    [
        e * (a * q[0] * q[0] - 2.0 * q[0] * a.powi(3) / (3.0 * r) + a.powi(5) / (5.0 * r * r)),
        4.0 * source.reduced_shear_modulus_pa() * a * (q[1] * q[1] + q[2] * q[2]),
        (8.0 / 3.0) * source.torsional_reduced_shear_modulus_pa() * a.powi(3) * q[3] * q[3],
        -PI * w * a * a,
    ]
}

#[test]
fn manufactured_joint_fold_stops_both_directions_and_releases_all_original_storage_once() {
    let source = source();
    let gate = CancelGate::new_clock_free();
    let (a, direction) = directions(&source, 4.0, 8.0);
    let delta = -a * a / source.normal_source().effective_radius_m();
    let loading = driver(&source, delta, direction, &gate);
    let limit = loading.limit().unwrap();
    assert_eq!(limit.kind, JkrCombinedBoundary::RadialFold);
    close(limit.contact_radius_m, a, 5e-13);
    close(limit.load_factor, 1.0, 5e-13);
    assert!(limit.radius_bracket_m[0] <= limit.contact_radius_m);
    assert!(limit.radius_bracket_m[1] >= limit.contact_radius_m);
    assert!(limit.radial_residual_n.abs() <= limit.radial_residual_tolerance_n);
    assert!(limit.radial_curvature_n_per_m.abs() <= limit.radial_curvature_tolerance_n_per_m);
    assert!(loading.preparation_evaluations() <= controls().maximum_evaluations);
    let initial = loading.initial_state(0.0, &gate).unwrap();
    let mut releases = [0.0; 2];
    for (index, sign) in [-1.0, 1.0].into_iter().enumerate() {
        let stopped = loading
            .advance_to_load_factor(initial, sign * 2.0 * limit.load_factor, &gate)
            .unwrap();
        let fold = stopped.state().contact().unwrap();
        let report = stopped.report();
        assert_eq!(stopped.previous_state(), initial);
        assert_eq!(report.boundary, Some(limit));
        assert_eq!(report.accepted_load_factor, sign * limit.load_factor);
        assert_eq!(report.untraversed_load_factor, sign * limit.load_factor);
        assert_eq!(fold.boundary, Some(JkrCombinedBoundary::RadialFold));
        assert!(fold.coupled_stiffness.is_none());
        assert_eq!(fold.displacement_radius_curvature_n_per_m, 0.0);
        assert_eq!(report.released_energy_j, 0.0);
        assert!(report.prescribed_work_j > 0.0);
        let before = potential(&source, initial.contact().unwrap());
        let after = potential(&source, fold);
        for (actual, expected) in [
            report.normal_elastic_change_j,
            report.tangential_elastic_change_j,
            report.torsional_elastic_change_j,
            report.surface_energy_change_j,
        ]
        .into_iter()
        .zip(core::array::from_fn::<_, 4, _>(|i| after[i] - before[i]))
        {
            close(actual, expected, 4e-12);
        }
        assert!(after[1] > 0.0 && after[2] > 0.0);
        close(
            report.prescribed_work_j,
            after.iter().sum::<f64>() - before.iter().sum::<f64>(),
            4e-12,
        );
        assert!(report.energy_residual_j.abs() <= report.energy_tolerance_j);

        // Eliminate U_a=U_aa=0 independently from the full potential. This
        // includes shear/twist release and is not the pure-normal snap value.
        let normal = source.normal_source();
        let e = normal.normal_state().elastic().reduced_modulus_pa();
        let a = fold.contact_radius_m;
        let depth = a * a / normal.effective_radius_m();
        let expected = e * a * (6.0 * depth * depth / 5.0 - 2.0 * delta * depth / 3.0)
            + (8.0 / 3.0)
                * source.torsional_reduced_shear_modulus_pa()
                * a.powi(3)
                * fold.coordinates[3].powi(2);
        let detached = loading.detach_at_fold(stopped.state(), &gate).unwrap();
        releases[index] = detached.report().released_energy_j;
        close(releases[index], expected, 4e-12);
        assert_eq!(detached.report().prescribed_work_j, 0.0);
        assert_eq!(detached.report().energy_residual_j, 0.0);
        assert_eq!(detached.state().contact(), None);
        assert_eq!(
            detached.state().coordinates(),
            stopped.state().coordinates()
        );
        assert_eq!(detached.state().stored_energy_j(), 0.0);
        assert_eq!(detached.state().released_energy_j(), releases[index]);
        assert!(matches!(
            loading.detach_at_fold(detached.state(), &gate),
            Err(JkrCombinedLoadingError::NotAtSeparatingFold)
        ));
        let free = loading
            .advance_to_load_factor(detached.state(), -2.0 * sign * limit.load_factor, &gate)
            .unwrap();
        assert_eq!(free.state().contact(), None);
        assert_eq!(free.report().prescribed_work_j, 0.0);
        assert_eq!(free.report().released_energy_j, 0.0);
        assert_eq!(free.report().untraversed_load_factor, 0.0);
        assert_eq!(free.state().released_energy_j(), releases[index]);
        assert_eq!(free.state().steps(), initial.steps() + 3);
    }
    assert_eq!(releases[0], releases[1]);
}

#[test]
fn combined_prescribed_work_matches_force_torque_quadrature_and_reversal_returns_storage() {
    let source = source();
    let gate = CancelGate::new_clock_free();
    let (a, direction) = directions(&source, 4.0, 8.0);
    let delta = -a * a / source.normal_source().effective_radius_m();
    let loading = driver(&source, delta, direction, &gate);
    let limit = loading.limit().unwrap().load_factor;
    let start = loading.initial_state(0.15 * limit, &gate).unwrap();
    let forward = loading
        .advance_to_load_factor(start, 0.7 * limit, &gate)
        .unwrap();
    let increment = forward.state().load_factor() - start.load_factor();
    let intervals = 256;
    let mut integral = 0.0;
    for i in 0..=intervals {
        let load = start.load_factor() + increment * i as f64 / intervals as f64;
        let point = source
            .evaluate_combined_displacement_controlled([
                delta,
                load * direction[0],
                load * direction[1],
                load * direction[2],
            ])
            .unwrap();
        let conjugate = point.force[1] * direction[0]
            + point.force[2] * direction[1]
            + point.force[3] * direction[2];
        let weight = if i == 0 || i == intervals {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        integral += weight * conjugate;
        assert!(point.contact_radius_m >= forward.report().minimum_contact_radius_m);
        assert!(point.contact_radius_m <= forward.report().maximum_contact_radius_m);
    }
    integral *= increment / (3.0 * intervals as f64);
    close(forward.report().prescribed_work_j, integral, 3e-9);
    let reverse = loading
        .advance_to_load_factor(forward.state(), start.load_factor(), &gate)
        .unwrap();
    close(
        reverse.report().prescribed_work_j,
        -forward.report().prescribed_work_j,
        2e-14,
    );
    assert_eq!(reverse.state().contact(), start.contact());
    let crossing = loading
        .advance_to_load_factor(forward.state(), -0.4 * limit, &gate)
        .unwrap();
    assert_eq!(
        crossing.report().maximum_contact_radius_m,
        loading
            .initial_state(0.0, &gate)
            .unwrap()
            .contact()
            .unwrap()
            .contact_radius_m
            .next_up()
    );
    assert_eq!(crossing.report().released_energy_j, 0.0);
    // The cached zero endpoint proves the reversal without a new source query.
    assert_eq!(crossing.report().evaluations, 1);
    let tiny = loading
        .advance_to_load_factor(start, start.load_factor().next_up(), &gate)
        .unwrap();
    assert!(tiny.report().prescribed_work_j > 0.0);
    let point = start.contact().unwrap();
    // At an adjacent load-factor bit, multiplication rounds the three actual
    // SI coordinate increments separately. Compare their actual force work,
    // not an unrepresentable ideal direction times the scalar increment.
    let after = tiny.state().coordinates();
    let expected = (1..4)
        .map(|i| point.force[i] * (after[i] - point.coordinates[i]))
        .sum::<f64>();
    close(tiny.report().prescribed_work_j, expected, 1e-8);
}

#[test]
fn opening_competes_with_fold_and_pure_mode_limits_keep_original_meaning() {
    let source = source();
    let gate = CancelGate::new_clock_free();
    for (b, t, kind) in [
        (8.0, 2.0, JkrCombinedBoundary::OpeningEdge),
        (2.0, 8.0, JkrCombinedBoundary::RadialFold),
    ] {
        let (a, direction) = directions(&source, b, t);
        let delta = a * a / source.normal_source().effective_radius_m();
        let loading = driver(&source, delta, direction, &gate);
        let limit = loading.limit().unwrap();
        assert_eq!(limit.kind, kind);
        let start = loading.initial_state(0.0, &gate).unwrap();
        let stopped = loading
            .advance_to_load_factor(start, 2.0 * limit.load_factor, &gate)
            .unwrap();
        assert_eq!(stopped.state().contact().unwrap().boundary, Some(kind));
        assert_eq!(stopped.report().released_energy_j, 0.0);
        assert!(matches!(
            loading.detach_at_fold(stopped.state(), &gate),
            Err(JkrCombinedLoadingError::NotAtSeparatingFold)
        ));
        let returned = loading
            .advance_to_load_factor(stopped.state(), 0.0, &gate)
            .unwrap();
        assert_eq!(returned.state().contact(), start.contact());
        close(
            returned.report().prescribed_work_j,
            -stopped.report().prescribed_work_j,
            2e-14,
        );
        if kind == JkrCombinedBoundary::OpeningEdge {
            close(limit.contact_radius_m, a, 5e-13);
            close(limit.load_factor, (16.0_f64 / (b + t)).sqrt(), 5e-13);
        }
    }
    let (a, direction) = directions(&source, 4.0, 8.0);
    let delta = -a * a / source.normal_source().effective_radius_m();
    let shear = driver(&source, delta, [direction[0], direction[1], 0.0], &gate);
    close(
        shear.limit().unwrap().load_factor * direction[0].hypot(direction[1]),
        source.critical_tangential_displacement_m(delta).unwrap(),
        5e-13,
    );
    let twist = driver(&source, delta, [0.0, 0.0, direction[2]], &gate);
    let original = source.torsional_displacement_limit(delta).unwrap().unwrap();
    assert_eq!(original.kind, JkrTorsionalLimitKind::RadialFold);
    close(
        twist.limit().unwrap().load_factor * direction[2],
        original.twist_angle_rad,
        5e-13,
    );
    close(
        twist.limit().unwrap().contact_radius_m,
        original.contact_radius_m,
        5e-13,
    );
    let zero = driver(&source, 0.0, [0.0, 0.0, direction[2]], &gate);
    assert_eq!(zero.limit(), None);
    let from = zero.initial_state(0.0, &gate).unwrap();
    let moved = zero.advance_to_load_factor(from, 0.5, &gate).unwrap();
    assert_eq!(moved.report().boundary, None);
    assert!(
        moved.state().contact().unwrap().contact_radius_m
            < from.contact().unwrap().contact_radius_m
    );
    assert!(driver(&source, 0.0, direction, &gate).limit().is_some());
}

#[test]
fn source_envelope_budget_identity_and_cancel_refusals_leave_exact_retry() {
    let source = source();
    let gate = CancelGate::new_clock_free();
    let (a, direction) = directions(&source, 4.0, 8.0);
    let delta = -a * a / source.normal_source().effective_radius_m();
    let loading = driver(&source, delta, direction, &gate);
    let start = loading.initial_state(0.2, &gate).unwrap();
    let expected = loading.advance_to_load_factor(start, -0.4, &gate).unwrap();
    for policy in [
        JkrCombinedLoadingControls {
            maximum_evaluations: 1,
            ..controls()
        },
        JkrCombinedLoadingControls {
            maximum_iterations: 1,
            ..controls()
        },
    ] {
        assert!(matches!(
            JkrCombinedLoading::new(
                &source,
                delta,
                [direction[0], direction[1]],
                direction[2],
                policy,
                &gate
            ),
            Err(JkrCombinedLoadingError::Budget)
        ));
    }
    let changed = driver(
        &source,
        delta,
        [direction[0], -direction[1], direction[2]],
        &gate,
    );
    assert!(matches!(
        changed.advance_to_load_factor(start, 0.3, &gate),
        Err(JkrCombinedLoadingError::StateMismatch)
    ));
    let policy = JkrCombinedLoading::new(
        &source,
        delta,
        [direction[0], direction[1]],
        direction[2],
        JkrCombinedLoadingControls {
            maximum_evaluations: 193,
            ..controls()
        },
        &gate,
    )
    .unwrap();
    assert!(matches!(
        policy.advance_to_load_factor(start, 0.3, &gate),
        Err(JkrCombinedLoadingError::StateMismatch)
    ));
    let restricted = SourceJkrMixedMode::new(
        source.normal_source().clone(),
        JkrMixedModeEnvelope {
            maximum_tangential_displacement_to_radius: 0.003,
            ..source.envelope()
        },
    )
    .unwrap();
    let limited = driver(&restricted, delta, direction, &gate);
    let limited_start = limited.initial_state(0.0, &gate).unwrap();
    assert!(
        limited
            .advance_to_load_factor(limited_start, 0.05, &gate)
            .is_ok()
    );
    assert!(matches!(
        limited.advance_to_load_factor(limited_start, 2.0, &gate),
        Err(JkrCombinedLoadingError::Source(_))
    ));
    assert_eq!(limited.initial_state(0.0, &gate).unwrap(), limited_start);
    assert!(matches!(
        loading.detach_at_fold(start, &gate),
        Err(JkrCombinedLoadingError::NotAtSeparatingFold)
    ));
    assert!(matches!(
        loading.initial_state(2.0 * loading.limit().unwrap().load_factor, &gate),
        Err(JkrCombinedLoadingError::BeyondFirstBoundary)
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        loading.advance_to_load_factor(start, 0.4, &cancelled),
        Err(JkrCombinedLoadingError::Cancelled)
    ));
    for bad in [f64::NAN, f64::INFINITY] {
        assert!(loading.advance_to_load_factor(start, bad, &gate).is_err());
    }
    assert_eq!(
        loading.advance_to_load_factor(start, -0.4, &gate).unwrap(),
        expected
    );
    assert_eq!(loading.initial_state(0.2, &gate).unwrap(), start);
}

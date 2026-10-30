//! Actual 3-D bridge shape, first variations and finite-rate wetting history.
#[path = "axisymmetric_contact_line/common.rs"]
mod common;
use common::*;
use core::f64::consts::PI;
use fs_contact::capillary_bridge::axisymmetric::contact_line::*;
use fs_exec::CancelGate;

#[test]
fn g1_young_equilibrium_reduces_to_original_axisymmetric_owner() {
    let source = source(1.1, 2.0, 1.08, 1.12);
    let gate = CancelGate::new_clock_free();
    for q in [0.0, 0.4, 2.0, 12.0] {
        let equilibrium = source
            .original()
            .solve_at_suction(q, numerical(), &gate)
            .unwrap();
        let state = source.initialize_at_suction(q, 1.1, &gate).unwrap();
        let p = source.evaluate(state).unwrap();
        let e = equilibrium.profile();
        close(p.gap_m, e.gap_m, 2e-16);
        close(p.neck_radius_m, e.neck_radius_m, 2e-16);
        close(p.wetted_radius_m, e.contact_radii_m[0], 2e-16);
        close(p.normal_force_n, e.normal_force_n, 2e-15);
        close(p.surface_free_energy_j, e.surface_free_energy_j, 2e-18);
        close(p.footprint_energy_gradient_n, 0.0, 2e-18);
        let step = source.advance(state, state.gap_m(), 1e-4, &gate).unwrap();
        assert_eq!(step.regime, ContactLineRegime::Pinned);
        assert_eq!(step.state(), state);
        assert_eq!(step.contact_line_dissipation_j, 0.0);
        assert_eq!(step.work_on_liquid_j, 0.0);
        if q > 0.0 {
            let relaxed = source.initialize(p.gap_m, 1.1, &gate).unwrap();
            close(relaxed.wetted_radius_m(), p.wetted_radius_m, 2e-16);
        }
    }
}

/// Independent RK4 of the radial-graph Young--Laplace ODE in z/a. No
/// production first-integral radius, angle quadrature or inverse solve is used.
fn original_ode(q: f64, height: f64, panels: usize) -> [f64; 4] {
    let derivative = |y: [f64; 4]| {
        let metric = 1.0 + y[1] * y[1];
        [
            y[1],
            metric / y[0] + q * metric * metric.sqrt(),
            y[0] * y[0],
            y[0] * metric.sqrt(),
        ]
    };
    let add =
        |x: [f64; 4], y: [f64; 4], factor: f64| core::array::from_fn(|i| x[i] + factor * y[i]);
    let ds = height / panels as f64;
    let mut y = [1.0, 0.0, 0.0, 0.0];
    for _ in 0..panels {
        let k1 = derivative(y);
        let k2 = derivative(add(y, k1, 0.5 * ds));
        let k3 = derivative(add(y, k2, 0.5 * ds));
        let k4 = derivative(add(y, k3, ds));
        y = core::array::from_fn(|i| y[i] + ds / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]));
    }
    y
}

#[test]
fn g1_actual_young_laplace_volume_area_pressure_and_first_variations() {
    let source = source(1.1, 2.0, 1.08, 1.12);
    let gate = CancelGate::new_clock_free();
    for q in [0.4, 2.0, 12.0] {
        let state = source.initialize_at_suction(q, 1.24, &gate).unwrap();
        let p = source.evaluate(state).unwrap();
        let a = p.neck_radius_m;
        let ode = original_ode(q, p.gap_m / (2.0 * a), 2048);
        close(a * ode[0], p.wetted_radius_m, 3e-14);
        close(2.0 * PI * a.powi(3) * ode[2], VOLUME, 2e-20);
        close(4.0 * PI * a * a * ode[3], p.liquid_gas_area_m2, 2e-16);
        close(ode[1], 1.24_f64.cos() / 1.24_f64.sin(), 2e-11);
        let plate_force = p.liquid_minus_gas_pressure_pa * PI * p.wetted_radius_m.powi(2)
            - 2.0 * PI * SIGMA * p.wetted_radius_m * p.contact_angle_radians.sin();
        close(plate_force, p.normal_force_n, 2e-15);
        let wet = source.original().surface_laws()[0].wetted_surface_energy_difference_j_m2();
        let energy = SIGMA * 4.0 * PI * a * a * ode[3] + 2.0 * wet * PI * (a * ode[0]).powi(2);
        close(energy, p.surface_free_energy_j, 2e-17);
        for index in 0..=8 {
            let point = source
                .meridian_point(state, index as f64 / 8.0, &gate)
                .unwrap();
            close(
                SIGMA * (point.meridional_curvature_per_m + point.azimuthal_curvature_per_m),
                p.liquid_minus_gas_pressure_pa,
                2e-10,
            );
        }
        for coordinate in 0..2 {
            let x = [p.gap_m, p.wetted_radius_m];
            // At high suction, fixed-volume angle sensitivity grows strongly
            // with radius/gap. A small physical perturbation avoids measuring
            // its finite nonlinear excursion as a first-variation defect.
            let delta = 1e-6 * x[0].min(x[1]);
            let mut lo = x;
            let mut hi = x;
            lo[coordinate] -= delta;
            hi[coordinate] += delta;
            let lower = source.initialize_at_radius(lo[0], lo[1], &gate).unwrap();
            let upper = source.initialize_at_radius(hi[0], hi[1], &gate).unwrap();
            let dl = source.evaluate(lower).unwrap();
            let du = source.evaluate(upper).unwrap();
            close(
                (du.surface_free_energy_j - dl.surface_free_energy_j) / (2.0 * delta),
                if coordinate == 0 {
                    -p.normal_force_n
                } else {
                    p.footprint_energy_gradient_n
                },
                2e-10,
            );
            close(
                (du.normal_force_n - dl.normal_force_n) / (2.0 * delta),
                if coordinate == 0 {
                    p.normal_force_gap_derivative_n_per_m
                } else {
                    p.normal_force_radius_derivative_n_per_m
                },
                2e-5,
            );
        }
    }
}

#[test]
fn g1_circular_perimeter_friction_and_pinning_balance_real_surface_work() {
    let source = source(1.1, 2.0, 1.08, 1.12);
    let gate = CancelGate::new_clock_free();
    for (theta, regime) in [
        (0.99, ContactLineRegime::Receding),
        (1.24, ContactLineRegime::Advancing),
    ] {
        let initial = source
            .initialize_at_suction(INITIAL_Q, theta, &gate)
            .unwrap();
        let mut state = initial;
        let mut loss = 0.0;
        let mut total_work = 0.0;
        for index in 0..40 {
            let next_gap = state.gap_m() + 1e-8;
            let step = source.advance(state, next_gap, 2e-4, &gate).unwrap();
            if index == 0 {
                assert_eq!(step.regime, regime);
            }
            assert!(step.regime == regime || step.regime == ContactLineRegime::Pinned);
            let db = step.state().wetted_radius_m() - state.wetted_radius_m();
            let mean_radius = 0.5 * (step.state().wetted_radius_m() + state.wetted_radius_m());
            close(
                step.viscous_line_dissipation_j,
                4.0 * PI * 2.0 * mean_radius * db * db / 2e-4,
                2e-25,
            );
            let pinning_per_radius = if db > 0.0 {
                4.0 * PI * SIGMA * (1.1_f64.cos() - 1.12_f64.cos())
            } else {
                -4.0 * PI * SIGMA * (1.08_f64.cos() - 1.1_f64.cos())
            };
            close(
                step.pinning_dissipation_j,
                pinning_per_radius * mean_radius * db,
                2e-25,
            );
            if step.regime == ContactLineRegime::Pinned {
                assert_eq!(db, 0.0);
                assert_eq!(step.contact_line_dissipation_j, 0.0);
            } else {
                assert!(step.viscous_line_dissipation_j > 0.0 && step.pinning_dissipation_j > 0.0);
            }
            assert_eq!(step.profile.liquid_volume_m3, VOLUME);
            assert!(step.shape_evaluations <= config().maximum_shape_evaluations);
            assert!(
                step.quadrature_evaluations
                    <= step.shape_evaluations as u64
                        * (u64::from(config().shape.maximum_quadrature_panels) + 1)
            );
            loss += step.contact_line_dissipation_j;
            total_work += step.work_on_liquid_j;
            state = step.state();
            close(
                state.profile().surface_free_energy_j + loss,
                initial.profile().surface_free_energy_j + total_work,
                3e-16,
            );
        }
        assert!(loss > 0.0);
    }
    let pinned_source = common::source(1.1, 2.0, 1.0, 1.2);
    let pinned = pinned_source
        .initialize_at_suction(INITIAL_Q, 1.1, &gate)
        .unwrap();
    let step = pinned_source
        .advance(pinned, pinned.gap_m() + 2e-7, 1e-3, &gate)
        .unwrap();
    assert_eq!(step.regime, ContactLineRegime::Pinned);
    assert_eq!(step.state().wetted_radius_m(), pinned.wetted_radius_m());
    assert_eq!(step.contact_line_dissipation_j, 0.0);
    assert!(step.work_on_liquid_j > 0.0);
    assert_ne!(
        step.state().contact_angle_radians(),
        pinned.contact_angle_radians()
    );
}

#[test]
fn g3_implicit_mean_force_tangent_tracks_pinned_advancing_and_receding_motion() {
    let source = source(1.1, 2.0, 1.08, 1.12);
    let gate = CancelGate::new_clock_free();
    for theta in [1.0, 1.1, 1.24] {
        let state = source
            .initialize_at_suction(INITIAL_Q, theta, &gate)
            .unwrap();
        let h = state.gap_m() + 2e-7;
        let epsilon = 2e-9;
        let step = source.advance(state, h, 2e-4, &gate).unwrap();
        let lo = source.advance(state, h - epsilon, 2e-4, &gate).unwrap();
        let hi = source.advance(state, h + epsilon, 2e-4, &gate).unwrap();
        assert_eq!(step.regime, lo.regime);
        assert_eq!(step.regime, hi.regime);
        close(
            (hi.mean_normal_force_n - lo.mean_normal_force_n) / (2.0 * epsilon),
            step.mean_force_gap_derivative_n_per_m,
            3e-5,
        );
    }
}

/// Independent fixed-volume shooting in neck scale and pressure; the radial
/// graph ODE is integrated in z, so the production inverse-angle formulation
/// does not serve as its own refinement oracle.
fn shooting(h: f64, b: f64, neck_guess: f64) -> (f64, f64) {
    let residual = |x: [f64; 2]| {
        let a = neck_guess * x[0];
        let ode = original_ode(x[1], h / (2.0 * a), 256);
        (
            [
                a * ode[0] / b - 1.0,
                2.0 * PI * a.powi(3) * ode[2] / VOLUME - 1.0,
            ],
            ode,
        )
    };
    let mut x = [1.0, INITIAL_Q];
    for _ in 0..16 {
        let (f, ode) = residual(x);
        if f.iter().all(|v| v.abs() < 3e-13) {
            let a = neck_guess * x[0];
            let theta = 1.0_f64.atan2(ode[1]);
            let energy =
                SIGMA * 4.0 * PI * a * a * ode[3] - 2.0 * SIGMA * 1.1_f64.cos() * PI * b * b;
            return (theta, energy);
        }
        let mut jac = [[0.0; 2]; 2];
        for j in 0..2 {
            let delta = 1e-5 * x[j].abs().max(1.0);
            let mut lo = x;
            let mut hi = x;
            lo[j] -= delta;
            hi[j] += delta;
            let fl = residual(lo).0;
            let fh = residual(hi).0;
            for i in 0..2 {
                jac[i][j] = (fh[i] - fl[i]) / (2.0 * delta);
            }
        }
        let determinant = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
        x[0] -= (jac[1][1] * f[0] - jac[0][1] * f[1]) / determinant;
        x[1] -= (-jac[1][0] * f[0] + jac[0][0] * f[1]) / determinant;
    }
    panic!("independent Young--Laplace shooting did not converge")
}

#[test]
fn g3_time_refinement_matches_independent_radial_graph_and_line_motion_odes() {
    let source = source(1.1, 2.0, 1.08, 1.12);
    let gate = CancelGate::new_clock_free();
    let initial = source
        .initialize_at_suction(INITIAL_Q, 1.24, &gate)
        .unwrap();
    let h = initial.gap_m();
    let a = initial.profile().neck_radius_m;
    let time = 0.003;
    let speed = |b: f64| {
        let theta = shooting(h, b, a).0;
        assert!(theta > 1.12);
        SIGMA / 2.0 * (1.12_f64.cos() - theta.cos())
    };
    let dt = time / 128.0;
    let mut reference = initial.wetted_radius_m();
    for _ in 0..128 {
        let k1 = speed(reference);
        let k2 = speed(reference + 0.5 * dt * k1);
        let k3 = speed(reference + 0.5 * dt * k2);
        let k4 = speed(reference + dt * k3);
        reference += dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }
    let reference_loss = shooting(h, initial.wetted_radius_m(), a).1 - shooting(h, reference, a).1;
    let mut errors = [0.0; 3];
    for (index, count) in [6, 12, 24].into_iter().enumerate() {
        let mut state = initial;
        let mut loss = 0.0;
        for _ in 0..count {
            let step = source
                .advance(state, h, time / count as f64, &gate)
                .unwrap();
            loss += step.contact_line_dissipation_j;
            state = step.state();
        }
        errors[index] = (state.wetted_radius_m() - reference).abs();
        if index == 2 {
            close(loss, reference_loss, 2e-13);
        }
    }
    assert!(
        errors[0] / errors[1] > 3.7 && errors[1] / errors[2] > 3.7,
        "{errors:?}"
    );
    assert!(errors[2] < 1e-9, "{errors:?}");
}

#[test]
fn g0_g4_symmetry_source_domain_accuracy_work_budget_and_cancellation_are_explicit() {
    let gate = CancelGate::new_clock_free();
    assert!(
        SourceAxisymmetricContactLineBridge::new(original(1.1, [2.0, 3.0], 1.08, 1.12), config())
            .is_err()
    );
    let source = source(1.1, 2.0, 1.08, 1.12);
    let state = source
        .initialize_at_suction(INITIAL_Q, 1.24, &gate)
        .unwrap();
    let changed = common::source(1.1, 3.0, 1.08, 1.12);
    assert!(matches!(
        changed.evaluate(state),
        Err(ContactLineError::StateMismatch)
    ));
    let expected = source.advance(state, state.gap_m(), 1e-4, &gate).unwrap();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        source.advance(state, state.gap_m(), 1e-4, &cancelled),
        Err(ContactLineError::Cancelled)
    ));
    assert!(
        source
            .advance(state, state.gap_m() + 1e-3, 1e-4, &gate)
            .is_err()
    );
    assert!(source.initialize_at_suction(-1.0, 1.1, &gate).is_err());
    assert!(source.initialize_at_suction(2.0, 1.55, &gate).is_err());
    assert!(source.advance(state, state.gap_m(), 0.1, &gate).is_err());
    assert_eq!(
        source.advance(state, state.gap_m(), 1e-4, &gate).unwrap(),
        expected
    );
    let mut small = config();
    small.maximum_shape_evaluations = 1;
    let limited =
        SourceAxisymmetricContactLineBridge::new(original(1.1, [2.0; 2], 1.08, 1.12), small)
            .unwrap();
    let old = limited
        .initialize_at_suction(INITIAL_Q, 1.24, &gate)
        .unwrap();
    assert!(matches!(
        limited.advance(old, old.gap_m(), 1e-4, &gate),
        Err(ContactLineError::Budget)
    ));
    assert_eq!(limited.evaluate(old).unwrap(), *old.profile());
    let mut strict = config();
    strict.energy_absolute_tolerance_j = 1e-30;
    strict.energy_relative_tolerance = 1e-20;
    let unresolved =
        SourceAxisymmetricContactLineBridge::new(original(1.1, [2.0; 2], 1.08, 1.12), strict)
            .unwrap();
    let old = unresolved
        .initialize_at_suction(INITIAL_Q, 1.1, &gate)
        .unwrap();
    assert!(matches!(
        unresolved.advance(old, old.gap_m(), 1e-4, &gate),
        Err(ContactLineError::Energy { .. })
    ));
    let mut strict = config();
    strict.force_absolute_tolerance_n = 1e-30;
    strict.force_relative_tolerance = 1e-20;
    let unresolved =
        SourceAxisymmetricContactLineBridge::new(original(1.1, [2.0; 2], 1.08, 1.12), strict)
            .unwrap();
    let old = unresolved
        .initialize_at_suction(INITIAL_Q, 1.1, &gate)
        .unwrap();
    assert!(matches!(
        unresolved.advance(old, old.gap_m(), 1e-4, &gate),
        Err(ContactLineError::Accuracy { .. })
    ));
}

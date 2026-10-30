//! Unlike circular footprints on the original three-dimensional liquid bridge.
#[path = "asymmetric_circular_contact_line/common.rs"]
mod common;
use common::*;
use core::f64::consts::PI;
use fs_contact::capillary_bridge::axisymmetric::asymmetric_contact_line::*;
use fs_exec::CancelGate;

#[test]
fn g1_unlike_young_equilibrium_matches_both_original_meridian_halves() {
    let source = source([1.02, 1.2], [2.0, 5.0], [1.0, 1.17], [1.05, 1.23]);
    let gate = CancelGate::new_clock_free();
    assert_ne!(
        source.original().wetting()[0].resolved().identity(),
        source.original().wetting()[1].resolved().identity()
    );
    for q in [0.0, 0.4, 2.0, 12.0] {
        let original = source
            .original()
            .solve_at_suction(q, numerical(), &gate)
            .unwrap();
        let state = source.initialize_at_suction(q, [1.02, 1.2], &gate).unwrap();
        let p = source.evaluate(state).unwrap();
        let e = original.profile();
        close(p.gap_m, e.gap_m, 2e-16);
        close(p.neck_radius_m, e.neck_radius_m, 2e-16);
        close(p.neck_height_m, e.neck_height_m, 2e-16);
        assert_ne!(p.wetted_radii_m[0], p.wetted_radii_m[1]);
        for i in 0..2 {
            close(p.wetted_radii_m[i], e.contact_radii_m[i], 2e-16);
            close(p.footprint_energy_gradients_n[i], 0.0, 2e-18);
        }
        close(p.normal_force_n, e.normal_force_n, 2e-15);
        close(p.surface_free_energy_j, e.surface_free_energy_j, 2e-18);
        let step = source.advance(state, state.gap_m(), 1e-4, &gate).unwrap();
        assert_eq!(step.regimes, [ContactLineRegime::Pinned; 2]);
        assert_eq!(step.state(), state);
        assert_eq!(step.contact_line_dissipation_j, 0.0);
        assert_eq!(step.work_on_liquid_j, 0.0);
    }
}

/// Independent radial-graph RK4: y=(r/a,dr/dz,int r² dz,int r sqrt(1+r'²) dz).
fn radial(q: f64, height: f64, panels: usize) -> [f64; 4] {
    let derivative = |y: [f64; 4]| {
        let m = 1.0 + y[1] * y[1];
        [
            y[1],
            m / y[0] + q * m * m.sqrt(),
            y[0] * y[0],
            y[0] * m.sqrt(),
        ]
    };
    let add = |x: [f64; 4], y: [f64; 4], s: f64| core::array::from_fn(|i| x[i] + s * y[i]);
    let dt = height / panels as f64;
    let mut y = [1.0, 0.0, 0.0, 0.0];
    for _ in 0..panels {
        let k1 = derivative(y);
        let k2 = derivative(add(y, k1, 0.5 * dt));
        let k3 = derivative(add(y, k2, 0.5 * dt));
        let k4 = derivative(add(y, k3, dt));
        y = core::array::from_fn(|i| y[i] + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]));
    }
    y
}

#[test]
fn g1_independent_young_laplace_geometry_and_complete_surface_hessian() {
    let source = source([1.05, 1.15], [2.0, 4.0], [1.0, 1.1], [1.1, 1.2]);
    let gate = CancelGate::new_clock_free();
    for q in [0.4, 2.0, 12.0] {
        let state = source
            .initialize_at_suction(q, [1.17, 1.28], &gate)
            .unwrap();
        let p = source.evaluate(state).unwrap();
        let a = p.neck_radius_m;
        let lower = radial(q, p.neck_height_m / a, 2048);
        let upper = radial(q, (p.gap_m - p.neck_height_m) / a, 2048);
        for (i, y) in [lower, upper].into_iter().enumerate() {
            close(a * y[0], p.wetted_radii_m[i], 3e-14);
            close(1.0_f64.atan2(y[1]), p.contact_angles_radians[i], 2e-11);
            let plate_force = p.liquid_minus_gas_pressure_pa * PI * p.wetted_radii_m[i].powi(2)
                - 2.0 * PI * SIGMA * p.wetted_radii_m[i] * p.contact_angles_radians[i].sin();
            close(plate_force, p.normal_force_n, 2e-15);
        }
        close(PI * a.powi(3) * (lower[2] + upper[2]), VOLUME, 2e-20);
        close(
            2.0 * PI * a * a * (lower[3] + upper[3]),
            p.liquid_gas_area_m2,
            2e-16,
        );
        let mut energy = SIGMA * 2.0 * PI * a * a * (lower[3] + upper[3]);
        for i in 0..2 {
            energy += source.original().surface_laws()[i].wetted_surface_energy_difference_j_m2()
                * PI
                * p.wetted_radii_m[i].powi(2);
        }
        close(energy, p.surface_free_energy_j, 2e-17);
        for fraction in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let point = source.meridian_point(state, fraction, &gate).unwrap();
            close(
                SIGMA * (point.meridional_curvature_per_m + point.azimuthal_curvature_per_m),
                p.liquid_minus_gas_pressure_pa,
                2e-10,
            );
        }
        let x = [p.gap_m, p.wetted_radii_m[0], p.wetted_radii_m[1]];
        let gradients = [
            -p.normal_force_n,
            p.footprint_energy_gradients_n[0],
            p.footprint_energy_gradients_n[1],
        ];
        for coordinate in 0..3 {
            let delta = 1e-6 * x[0].min(x[1]).min(x[2]);
            let mut lo = x;
            let mut hi = x;
            lo[coordinate] -= delta;
            hi[coordinate] += delta;
            let low = source
                .initialize_at_radii(lo[0], [lo[1], lo[2]], &gate)
                .unwrap();
            let high = source
                .initialize_at_radii(hi[0], [hi[1], hi[2]], &gate)
                .unwrap();
            let l = source.evaluate(low).unwrap();
            let h = source.evaluate(high).unwrap();
            close(
                (h.surface_free_energy_j - l.surface_free_energy_j) / (2.0 * delta),
                gradients[coordinate],
                3e-10,
            );
            let gl = [
                -l.normal_force_n,
                l.footprint_energy_gradients_n[0],
                l.footprint_energy_gradients_n[1],
            ];
            let gh = [
                -h.normal_force_n,
                h.footprint_energy_gradients_n[0],
                h.footprint_energy_gradients_n[1],
            ];
            for i in 0..3 {
                close(
                    (gh[i] - gl[i]) / (2.0 * delta),
                    p.energy_hessian_n_per_m[i][coordinate],
                    5e-5,
                );
            }
        }
    }
}

#[test]
fn g1_independent_pinning_and_actual_circular_perimeter_loss_close_work() {
    let theta = [1.05, 1.15];
    let zeta = [2.0, 4.0];
    let rec = [1.0, 1.1];
    let adv = [1.1, 1.2];
    let source = source(theta, zeta, rec, adv);
    let gate = CancelGate::new_clock_free();
    for (angles, expected) in [
        (
            [1.05, 1.28],
            [ContactLineRegime::Pinned, ContactLineRegime::Advancing],
        ),
        (
            [0.94, 1.15],
            [ContactLineRegime::Receding, ContactLineRegime::Pinned],
        ),
    ] {
        let mut state = source
            .initialize_at_suction(INITIAL_Q, angles, &gate)
            .unwrap();
        let initial = state;
        let mut total_loss = 0.0;
        for _ in 0..12 {
            let step = source.advance(state, state.gap_m(), 2e-5, &gate).unwrap();
            assert_eq!(step.regimes, expected);
            for i in 0..2 {
                let db = step.state().wetted_radii_m()[i] - state.wetted_radii_m()[i];
                let mean = 0.5 * (step.state().wetted_radii_m()[i] + state.wetted_radii_m()[i]);
                close(
                    step.per_line_viscous_dissipation_j[i],
                    2.0 * PI * zeta[i] * mean * db * db / 2e-5,
                    2e-24,
                );
                let threshold = if db > 0.0 {
                    2.0 * PI * SIGMA * (theta[i].cos() - adv[i].cos())
                } else if db < 0.0 {
                    -2.0 * PI * SIGMA * (rec[i].cos() - theta[i].cos())
                } else {
                    0.0
                };
                close(
                    step.per_line_pinning_dissipation_j[i],
                    threshold * mean * db,
                    2e-24,
                );
                if expected[i] == ContactLineRegime::Pinned {
                    assert_eq!(db, 0.0);
                }
            }
            assert!(step.energy_balance_residual_j.abs() + step.estimated_energy_error_j < 3e-17);
            total_loss += step.contact_line_dissipation_j;
            state = step.state();
        }
        close(
            initial.profile().surface_free_energy_j - state.profile().surface_free_energy_j,
            total_loss,
            4e-17,
        );
    }
}

#[test]
fn g3_exchanging_original_surfaces_exchanges_footprints_and_each_line_friction_matters() {
    let gate = CancelGate::new_clock_free();
    let source = source([1.05, 1.15], [2.0, 4.0], [1.0, 1.1], [1.1, 1.2]);
    let swapped = common::source([1.15, 1.05], [4.0, 2.0], [1.1, 1.0], [1.2, 1.1]);
    let s = source
        .initialize_at_suction(INITIAL_Q, [1.22, 1.28], &gate)
        .unwrap();
    let t = swapped
        .initialize_at_suction(INITIAL_Q, [1.28, 1.22], &gate)
        .unwrap();
    let a = source.advance(s, s.gap_m() + 1e-7, 2e-4, &gate).unwrap();
    let b = swapped.advance(t, t.gap_m() + 1e-7, 2e-4, &gate).unwrap();
    for i in 0..2 {
        close(
            a.state().wetted_radii_m()[i],
            b.state().wetted_radii_m()[1 - i],
            2e-15,
        );
        close(
            a.profile.contact_angles_radians[i],
            b.profile.contact_angles_radians[1 - i],
            2e-11,
        );
        close(
            a.per_line_viscous_dissipation_j[i],
            b.per_line_viscous_dissipation_j[1 - i],
            2e-20,
        );
    }
    close(
        a.profile.neck_height_m,
        b.profile.gap_m - b.profile.neck_height_m,
        2e-15,
    );
    close(a.mean_normal_force_n, b.mean_normal_force_n, 2e-14);
    close(
        a.contact_line_dissipation_j,
        b.contact_line_dissipation_j,
        2e-19,
    );
    let slow = common::source([1.05, 1.15], [2.0, 16.0], [1.0, 1.1], [1.1, 1.2]);
    let old = slow
        .initialize_at_suction(INITIAL_Q, [1.22, 1.28], &gate)
        .unwrap();
    let slower = slow.advance(old, old.gap_m() + 1e-7, 2e-4, &gate).unwrap();
    assert!(
        slower.state().wetted_radii_m()[1] - old.wetted_radii_m()[1]
            < 0.5 * (a.state().wetted_radii_m()[1] - s.wetted_radii_m()[1])
    );
    assert_ne!(source.identity(), slow.identity());
}

#[test]
fn g3_implicit_normal_tangent_retains_both_moving_and_pinned_footprints() {
    let source = source([1.05, 1.15], [2.0, 4.0], [1.0, 1.1], [1.1, 1.2]);
    let gate = CancelGate::new_clock_free();
    for angles in [[1.05, 1.15], [1.05, 1.28], [0.94, 1.28], [1.22, 1.28]] {
        let state = source
            .initialize_at_suction(INITIAL_Q, angles, &gate)
            .unwrap();
        let h = state.gap_m() + 2e-7;
        let epsilon = 2e-9;
        let step = source.advance(state, h, 2e-4, &gate).unwrap();
        let lo = source.advance(state, h - epsilon, 2e-4, &gate).unwrap();
        let hi = source.advance(state, h + epsilon, 2e-4, &gate).unwrap();
        assert_eq!(step.regimes, lo.regimes);
        assert_eq!(step.regimes, hi.regimes);
        close(
            (hi.mean_normal_force_n - lo.mean_normal_force_n) / (2.0 * epsilon),
            step.mean_force_gap_derivative_n_per_m,
            4e-5,
        );
    }
}

/// Independently solve for neck scale, pressure and neck position by radial
/// shooting. The numerical reference never calls the production angle inverse.
fn shooting(h: f64, b: [f64; 2], neck_guess: f64, neck_fraction: f64) -> ([f64; 2], f64) {
    let residual = |x: [f64; 3]| {
        let a = neck_guess * x[0];
        let lower = radial(x[1], h * x[2] / a, 256);
        let upper = radial(x[1], h * (1.0 - x[2]) / a, 256);
        (
            [
                a * lower[0] / b[0] - 1.0,
                a * upper[0] / b[1] - 1.0,
                PI * a.powi(3) * (lower[2] + upper[2]) / VOLUME - 1.0,
            ],
            [lower, upper],
        )
    };
    let mut x = [1.0, INITIAL_Q, neck_fraction];
    for _ in 0..20 {
        let (f, ys) = residual(x);
        if f.iter().all(|v| v.abs() < 3e-13) {
            let a = neck_guess * x[0];
            let angles = ys.map(|y| 1.0_f64.atan2(y[1]));
            let energy = SIGMA * 2.0 * PI * a * a * (ys[0][3] + ys[1][3])
                - SIGMA * PI * (1.05_f64.cos() * b[0] * b[0] + 1.15_f64.cos() * b[1] * b[1]);
            return (angles, energy);
        }
        let mut augmented = [[0.0; 4]; 3];
        for i in 0..3 {
            augmented[i][3] = f[i];
        }
        for j in 0..3 {
            let delta = 1e-5 * x[j].abs().max(1.0);
            let mut lo = x;
            let mut hi = x;
            lo[j] -= delta;
            hi[j] += delta;
            let fl = residual(lo).0;
            let fh = residual(hi).0;
            for i in 0..3 {
                augmented[i][j] = (fh[i] - fl[i]) / (2.0 * delta);
            }
        }
        for j in 0..3 {
            let row = (j..3)
                .max_by(|a, b| augmented[*a][j].abs().total_cmp(&augmented[*b][j].abs()))
                .unwrap();
            augmented.swap(row, j);
            let pivot = augmented[j][j];
            for k in j..4 {
                augmented[j][k] /= pivot;
            }
            for i in 0..3 {
                if i != j {
                    let factor = augmented[i][j];
                    for k in j..4 {
                        augmented[i][k] -= factor * augmented[j][k];
                    }
                }
            }
        }
        for i in 0..3 {
            x[i] -= augmented[i][3];
        }
    }
    panic!("independent two-ended radial shooting did not converge")
}

#[test]
fn g3_second_order_refinement_against_independent_shape_and_two_line_motion_odes() {
    let source = source([1.05, 1.15], [2.0, 4.0], [1.0, 1.1], [1.1, 1.2]);
    let gate = CancelGate::new_clock_free();
    let initial = source
        .initialize_at_suction(INITIAL_Q, [1.22, 1.28], &gate)
        .unwrap();
    let h = initial.gap_m();
    let a = initial.profile().neck_radius_m;
    let fraction = initial.profile().neck_height_m / h;
    let speed = |b: [f64; 2]| {
        let angles = shooting(h, b, a, fraction).0;
        assert!(angles[0] > 1.1 && angles[1] > 1.2, "{angles:?}");
        [
            SIGMA / 2.0 * (1.1_f64.cos() - angles[0].cos()),
            SIGMA / 4.0 * (1.2_f64.cos() - angles[1].cos()),
        ]
    };
    let time = 0.0015;
    let dt = time / 128.0;
    let add = |x: [f64; 2], y: [f64; 2], s: f64| core::array::from_fn(|i| x[i] + s * y[i]);
    let mut reference = initial.wetted_radii_m();
    for _ in 0..128 {
        let k1 = speed(reference);
        let k2 = speed(add(reference, k1, 0.5 * dt));
        let k3 = speed(add(reference, k2, 0.5 * dt));
        let k4 = speed(add(reference, k3, dt));
        reference = core::array::from_fn(|i| {
            reference[i] + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])
        });
    }
    let reference_loss = shooting(h, initial.wetted_radii_m(), a, fraction).1
        - shooting(h, reference, a, fraction).1;
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
        errors[index] = (state.wetted_radii_m()[0] - reference[0])
            .hypot(state.wetted_radii_m()[1] - reference[1]);
        if index == 2 {
            close(loss, reference_loss, 2e-13);
        }
    }
    assert!(
        errors[0] / errors[1] > 3.6 && errors[1] / errors[2] > 3.6,
        "{errors:?}"
    );
    assert!(errors[2] < 1e-9, "{errors:?}");
}

#[test]
fn g0_g4_original_source_accuracy_domain_budget_and_cancellation_remain_explicit() {
    use fs_contact::capillary_bridge::{BridgeError, axisymmetric::AxisymmetricError};
    let gate = CancelGate::new_clock_free();
    let source = source([1.05, 1.15], [2.0, 4.0], [1.0, 1.1], [1.1, 1.2]);
    let old = source
        .initialize_at_suction(INITIAL_Q, [1.22, 1.28], &gate)
        .unwrap();
    let expected = source.advance(old, old.gap_m(), 1e-4, &gate).unwrap();
    assert_ne!(
        source.original().wetting()[0].context().history,
        source.original().wetting()[1].context().history
    );
    assert_eq!(
        source.original().wetting()[0].context().history,
        "original surface 0"
    );
    assert_eq!(
        source.original().wetting()[1].context().history,
        "original surface 1"
    );
    assert!(matches!(
        different_environment(),
        Err(AxisymmetricError::Physical(BridgeError::Invalid {
            field: "common original surface environment"
        }))
    ));
    let other = common::source([1.05, 1.15], [2.0, 5.0], [1.0, 1.1], [1.1, 1.2]);
    assert!(matches!(
        other.evaluate(old),
        Err(ContactLineError::StateMismatch)
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        source.advance(old, old.gap_m(), 1e-4, &cancelled),
        Err(ContactLineError::Cancelled)
    ));
    assert!(
        source
            .initialize_at_suction(-1.0, [1.1, 1.2], &gate)
            .is_err()
    );
    assert!(
        source
            .initialize_at_suction(2.0, [1.1, 1.55], &gate)
            .is_err()
    );
    assert!(
        source
            .advance(old, old.gap_m() + 1e-3, 1e-4, &gate)
            .is_err()
    );
    assert!(source.advance(old, old.gap_m(), 0.1, &gate).is_err());
    assert_eq!(
        source.advance(old, old.gap_m(), 1e-4, &gate).unwrap(),
        expected
    );
    let mut small = config();
    small.maximum_shape_evaluations = 2;
    let limited = SourceAsymmetricCircularContactLineBridge::new(
        original([1.05, 1.15], [2.0, 4.0], [1.0, 1.1], [1.1, 1.2]),
        small,
    )
    .unwrap();
    let old = limited
        .initialize_at_suction(INITIAL_Q, [1.22, 1.28], &gate)
        .unwrap();
    assert!(matches!(
        limited.advance(old, old.gap_m(), 1e-4, &gate),
        Err(ContactLineError::Budget)
    ));
    assert_eq!(limited.evaluate(old).unwrap(), *old.profile());
    for force in [false, true] {
        let mut strict = config();
        if force {
            strict.force_absolute_tolerance_n = 1e-30;
            strict.force_relative_tolerance = 1e-20;
        } else {
            strict.energy_absolute_tolerance_j = 1e-30;
            strict.energy_relative_tolerance = 1e-20;
        }
        let unresolved = SourceAsymmetricCircularContactLineBridge::new(
            original([1.05, 1.15], [2.0, 4.0], [1.0, 1.1], [1.1, 1.2]),
            strict,
        )
        .unwrap();
        let old = unresolved
            .initialize_at_suction(INITIAL_Q, [1.05, 1.15], &gate)
            .unwrap();
        let error = unresolved
            .advance(old, old.gap_m(), 1e-4, &gate)
            .unwrap_err();
        assert!(
            if force {
                matches!(
                    error,
                    ContactLineError::Accuracy { .. } | ContactLineError::NotConverged
                )
            } else {
                matches!(error, ContactLineError::Energy { .. })
            },
            "{error:?}"
        );
    }
}

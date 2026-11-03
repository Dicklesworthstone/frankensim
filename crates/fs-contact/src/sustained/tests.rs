use super::*;

fn scales() -> NaturalMapScales {
    NaturalMapScales {
        normal_n_per_m: 100.0,
        tangent_n_per_m: 50.0,
    }
}
fn point(force_n: [f64; 3], gap_m: f64, slip_increment_m: [f64; 2], mu: f64) -> ContactPoint {
    ContactPoint {
        force_n,
        gap_m,
        slip_increment_m,
        friction_coefficient: mu,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e} (tolerance {tolerance:e})"
    );
}
fn zero(value: [f64; 3]) {
    for v in value {
        close(v, 0.0, 1e-12);
    }
}

#[test]
fn normal_support_static_capacity_and_sliding_force_obey_the_original_inclusions() {
    // A 20 N normal reaction supplies a 5 N disk at mu=1/4.
    let stick = point([20.0, -4.0, 0.0], 0.0, [0.0; 2], 0.25);
    let slide = point([20.0, -5.0, 0.0], 0.0, [1.0, 0.0], 0.25);
    for p in [
        stick,
        slide,
        point([0.0; 3], 0.2, [0.3, -0.4], 0.25),
        point([20.0, 0.0, 0.0], 0.0, [1.0, -2.0], 0.0),
    ] {
        zero(linearize(p, scales()).unwrap().residual_n);
        // Numerical scaling does not alter the physical zero or create a spring.
        zero(
            linearize(
                p,
                NaturalMapScales {
                    normal_n_per_m: 1e6,
                    tangent_n_per_m: 1e3,
                },
            )
            .unwrap()
            .residual_n,
        );
    }
    close(-(slide.force_n[1] * slide.slip_increment_m[0]), 5.0, 1e-14);
    let too_large = linearize(point([20.0, -6.0, 0.0], 0.0, [0.0; 2], 0.25), scales()).unwrap();
    close(too_large.residual_n[1], -1.0, 1e-14);
    let wrong_direction =
        linearize(point([20.0, 5.0, 0.0], 0.0, [1.0, 0.0], 0.25), scales()).unwrap();
    close(wrong_direction.residual_n[1], 10.0, 1e-14);
    let separated_load = linearize(point([20.0, 0.0, 0.0], 0.1, [0.0; 2], 0.25), scales()).unwrap();
    close(separated_load.residual_n[0], 10.0, 1e-14);
    // Negative forces/gaps are valid nonlinear TRIALS, with nonzero residuals.
    assert!(
        linearize(point([-1.0, 0.0, 0.0], 0.0, [0.0; 2], 0.25), scales())
            .unwrap()
            .residual_n[0]
            < 0.0
    );
    assert!(
        linearize(point([20.0, 0.0, 0.0], -0.01, [0.0; 2], 0.25), scales())
            .unwrap()
            .residual_n[0]
            < 0.0
    );
}

#[test]
fn selected_kink_actions_are_fixed_linear_maps_for_native_krylov() {
    let df_a = [0.3, -0.4, 0.7];
    let df_b = [-0.2, 0.1, -0.6];
    let dp_a = [0.004, -0.003, 0.002];
    let dp_b = [-0.001, 0.005, -0.006];
    for p in [
        point([0.0; 3], 0.0, [0.0; 2], 0.25),
        point([20.0, 5.0, 0.0], 0.0, [0.0; 2], 0.25),
        point([20.0, -5.0, 0.0], 0.0, [1.0, 0.0], 0.25),
    ] {
        let fixed = linearize(p, scales()).unwrap();
        let a = fixed.action(df_a, dp_a).unwrap();
        let b = fixed.action(df_b, dp_b).unwrap();
        let summed = fixed
            .action(
                core::array::from_fn(|i| df_a[i] + df_b[i]),
                core::array::from_fn(|i| dp_a[i] + dp_b[i]),
            )
            .unwrap();
        let reversed = fixed.action(df_a.map(|x| -x), dp_a.map(|x| -x)).unwrap();
        for i in 0..3 {
            close(summed[i], a[i] + b[i], 1e-12);
            close(reversed[i], -a[i], 1e-12);
        }
        zero(fixed.action([0.0; 3], [0.0; 3]).unwrap());
    }
    let apex = linearize(point([0.0; 3], 0.0, [0.0; 2], 0.25), scales()).unwrap();
    assert_eq!(
        apex.force_derivative,
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    );
    assert_eq!(apex.path_derivative, [[0.0; 3]; 3]);
    let boundary = linearize(point([20.0, 5.0, 0.0], 0.0, [0.0; 2], 0.25), scales()).unwrap();
    close(boundary.force_derivative[1][0], -0.25, 1e-14);
    close(boundary.force_derivative[1][1], 1.0, 1e-14);
    close(boundary.path_derivative[2][2], 50.0, 1e-14);
}

#[test]
fn analytic_contact_blocks_match_independent_differences_away_from_mode_kinks() {
    let base = point([10.0, 2.0, -1.0], 0.001, [0.03, 0.02], 0.4);
    let scales = NaturalMapScales {
        normal_n_per_m: 1000.0,
        tangent_n_per_m: 500.0,
    };
    let df = [0.3, -0.7, 1.2];
    let dp = [-0.002, 0.004, -0.003];
    let h = 1e-6;
    let shift = |sign: f64| ContactPoint {
        force_n: core::array::from_fn(|i| base.force_n[i] + sign * h * df[i]),
        gap_m: base.gap_m + sign * h * dp[0],
        slip_increment_m: core::array::from_fn(|i| base.slip_increment_m[i] + sign * h * dp[i + 1]),
        ..base
    };
    let plus = linearize(shift(1.0), scales).unwrap().residual_n;
    let minus = linearize(shift(-1.0), scales).unwrap().residual_n;
    let expected = linearize(base, scales).unwrap().action(df, dp).unwrap();
    for i in 0..3 {
        close((plus[i] - minus[i]) / (2.0 * h), expected[i], 1e-7);
    }
    let slide = linearize(
        point([20.0, -5.0, 0.0], 0.0, [1.0, 0.0], 0.25),
        super::tests::scales(),
    )
    .unwrap();
    close(slide.force_derivative[1][0], 0.25, 1e-14);
    close(slide.force_derivative[2][2], 10.0 / 11.0, 1e-14);
    close(slide.path_derivative[2][2], 50.0 / 11.0, 1e-14);
}

#[test]
fn invalid_scales_and_unrepresentable_maps_cannot_become_contact_acceptance() {
    let p = point([20.0, -4.0, 0.0], 0.0, [0.0; 2], 0.25);
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            linearize(
                p,
                NaturalMapScales {
                    normal_n_per_m: invalid,
                    ..scales()
                }
            )
            .is_err()
        );
        assert!(
            linearize(
                p,
                NaturalMapScales {
                    tangent_n_per_m: invalid,
                    ..scales()
                }
            )
            .is_err()
        );
    }
    assert!(
        linearize(
            ContactPoint {
                friction_coefficient: -0.25,
                ..p
            },
            scales()
        )
        .is_err()
    );
    assert!(
        linearize(
            ContactPoint {
                gap_m: f64::NAN,
                ..p
            },
            scales()
        )
        .is_err()
    );
    assert!(
        linearize(
            ContactPoint {
                force_n: [f64::MAX, 0.0, 0.0],
                friction_coefficient: 2.0,
                ..p
            },
            scales()
        )
        .is_err()
    );
    assert!(
        linearize(
            ContactPoint {
                slip_increment_m: [f64::MAX; 2],
                ..p
            },
            scales()
        )
        .is_err()
    );
    let fixed = linearize(p, scales()).unwrap();
    assert!(fixed.action([f64::NAN, 0.0, 0.0], [0.0; 3]).is_err());
    assert!(fixed.action([0.0; 3], [f64::MAX; 3]).is_err());
}

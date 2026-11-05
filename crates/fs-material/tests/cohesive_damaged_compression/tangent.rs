//! G0/G1 derivatives of the original path measure, with independent references.
use super::*;

fn pressure(old: CohesivePoint<'_>, end: [f64; 3]) -> f64 {
    old.trial_over_time_with_damage_weighted_compression(end, 0.01, 1e-8)
        .unwrap()
        .2
}

fn gradient(old: CohesivePoint<'_>, end: [f64; 3]) -> [f64; 3] {
    let original = old
        .trial_over_time_with_damage_weighted_compression(end, 0.01, 1e-8)
        .unwrap();
    let mean_tangent = old
        .trial_over_time_with_mean_traction_tangent(end, 0.01, 1e-8)
        .unwrap();
    let result = old
        .trial_over_time_with_damage_weighted_compression_tangent(end, 0.01, 1e-8)
        .unwrap();
    assert_eq!((result.0, result.1, result.2), original);
    assert_eq!(result.3, mean_tangent.2);
    result.4.expect("smooth compression path gradient")
}

fn check_difference(old: CohesivePoint<'_>, end: [f64; 3]) {
    let saved = old;
    let exact = gradient(old, end);
    for j in 0..3 {
        let h = 1e-6 * (2e-4 + old.jump_m()[j].abs() + end[j].abs());
        let mut plus = end;
        let mut minus = end;
        plus[j] += h;
        minus[j] -= h;
        let difference = (pressure(old, plus) - pressure(old, minus)) / (2.0 * h);
        assert!(
            (exact[j] - difference).abs() <= 2.0 + 3e-6 * exact[j].abs().max(difference.abs()),
            "start={:?}, end={end:?}, column={j}: analytic={:e}, difference={difference:e}",
            old.jump_m(),
            exact[j],
        );
    }
    assert_eq!(old, saved);
}

#[test]
fn smooth_gradient_follows_original_onset_unloading_rotation_and_rupture() {
    for law in laws() {
        for (start, end) in [
            ([-2e-4, 0.0, 0.0], [-4e-4, 8e-4, 3e-4]),
            ([1.2e-4, 6e-5, -3e-5], [-4e-4, 0.0012, -5e-4]),
            ([1.2e-4, 6e-5, -3e-5], [-4e-4, 0.012, 0.004]),
            ([4e-4, 2e-4, 1e-4], [-4e-4, -0.003, 0.001]),
            ([-3e-4, 2e-4, 0.0], [4e-4, 0.012, -0.004]),
            ([-1e-4, 6e-4, -3e-4], [-3e-4, 1e-4, 2e-4]),
            ([-1e-4, 0.014, 0.0], [-4e-4, -0.003, 0.001]),
            ([-3e-4, 0.0, 0.0], [1e-4, 0.0, 0.0]),
            ([1e-4, 0.0, 0.0], [4e-4, 0.004, 0.0]),
            ([-2e-4, 5e-4, 1e-4], [0.0, 8e-4, -2e-4]),
            ([0.0, 5e-4, 1e-4], [-2e-4, 8e-4, -2e-4]),
        ] {
            let old = law.bonded_point().trial(start, 1e-8).unwrap().point;
            check_difference(old, end);
        }
    }
}

#[test]
fn first_bk_toughness_sensitivity_survives_reclosure_and_shear_basis_rotation() {
    let rotate = |g: [f64; 3]| [g[0], 0.6 * g[1] - 0.8 * g[2], 0.8 * g[1] + 0.6 * g[2]];
    for exponent in [1.0, 1.35, 2.0, 3.7] {
        let mut parameters = laws()[1].mixed_mode_parameters().unwrap();
        parameters.bk_exponent = exponent;
        let law = BilinearCohesive::new_mixed_mode(parameters, 2e9, 0.02).unwrap();
        let start = [1.2e-4, 6e-5, -3e-5];
        let old = law.bonded_point().trial(start, 1e-8).unwrap().point;
        let rotated = law.bonded_point().trial(rotate(start), 1e-8).unwrap().point;
        assert_eq!(old.damage(), 0.0);
        for end in [[-4e-4, 0.0012, -5e-4], [-4e-4, 0.012, 0.004]] {
            check_difference(old, end);
            let original_gradient = gradient(old, end);
            let rotated_gradient = gradient(rotated, rotate(end));
            for (actual, expected) in rotated_gradient.into_iter().zip(rotate(original_gradient)) {
                close(actual, expected, 5e-10);
            }
            close(pressure(old, end), pressure(rotated, rotate(end)), 2e-12);
        }
    }
}

#[test]
fn intact_and_fully_failed_limits_have_independent_closed_form_gradients() {
    for law in laws() {
        let intact = law
            .bonded_point()
            .trial([-2e-4, 0.0, 0.0], 1e-8)
            .unwrap()
            .point;
        let end = [-6e-4, 3e-5, -2e-5];
        assert_eq!(pressure(intact, end), 0.0);
        assert_eq!(gradient(intact, end), [0.0; 3]);

        let failed = law
            .bonded_point()
            .trial([-2e-4, 0.014, 0.0], 1e-8)
            .unwrap()
            .point;
        assert!(failed.failed());
        assert_eq!(gradient(failed, end), [-1e9, 0.0, 0.0]);
        close(pressure(failed, end), 2e9 * 4e-4, 1e-14);

        // A triangular compression interval ending in opening has an exact
        // moving-boundary derivative; no damage or cohesive area remains.
        let opening = [4e-4, -1e-4, 2e-5];
        let initial_compression: f64 = 2e-4;
        let length = initial_compression + opening[0];
        close(
            gradient(failed, opening)[0],
            -1e9 * (initial_compression / length).powi(2),
            2e-14,
        );
        assert_eq!(&gradient(failed, opening)[1..], &[0.0, 0.0]);
        close(
            pressure(failed, opening),
            1e9 * initial_compression.powi(2) / length,
            2e-14,
        );
    }
}

#[test]
fn newly_damaged_area_gradient_resolves_short_softening_without_subtracting_pressures() {
    let law = laws()[0].clone();
    let onset = law.onset_separation_m();
    let failure = law.failure_separation_m();
    let factor = failure / (failure - onset);
    let old = law
        .bonded_point()
        .trial([-2e-4, onset, 0.0], 1e-8)
        .unwrap()
        .point;
    for fraction in [1e-10, 1e-8, 1e-6, 1e-4] {
        let end = [-2e-4, onset * (1.0 + fraction), 0.0];
        let epsilon = (end[1] - onset) / onset;
        let actual = gradient(old, end);
        // For r(s)=onset*(1+epsilon*s), d(s)=factor*epsilon*s/(1+epsilon*s).
        // These independent scalar series integrate s*d(s) and the endpoint
        // shear derivative of d(s), retaining the tiny newly damaged area.
        let mut normal_moment = 0.0;
        let mut shear_moment = 0.0;
        for k in 0..8 {
            let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
            normal_moment += sign * epsilon.powi(k + 1) / f64::from(k + 3);
            shear_moment += sign * f64::from(k + 1) * epsilon.powi(k) / f64::from(k + 2);
        }
        close(actual[0], -2e9 * factor * normal_moment, 2e-6);
        close(
            actual[1],
            (2e9 * 2e-4 / onset) * factor * shear_moment,
            2e-11,
        );
        assert_eq!(actual[2], 0.0);
        assert!(pressure(old, end) > 0.0);
    }
}

#[test]
fn compression_gradient_has_no_bond_viscosity_or_step_length_dependence() {
    for law in laws() {
        let viscous = law
            .clone()
            .with_viscosity(CohesiveViscosity::new(1e5, 2e5, 100.0).unwrap());
        let start = [-2e-4, 0.0, 0.0];
        let end = [-4e-4, 8e-4, 3e-4];
        let old = law.bonded_point().trial(start, 1e-8).unwrap().point;
        let damped = viscous
            .bonded_point()
            .trial_over_time(start, 0.01, 1e-8)
            .unwrap()
            .equilibrium
            .point;
        let expected = gradient(old, end);
        let original = old
            .trial_over_time_with_damage_weighted_compression_tangent(end, 0.01, 1e-8)
            .unwrap();
        for dt in [0.005, 0.01, 0.04] {
            let response = damped
                .trial_over_time_with_damage_weighted_compression_tangent(end, dt, 1e-8)
                .unwrap();
            assert_eq!(response.2, original.2);
            assert_eq!(response.4, Some(expected));
            assert_ne!(response.3, original.3);
            assert!(response.0.viscous_dissipated_j_m2 > 0.0);
        }
    }
}

#[test]
fn unsupported_classical_gradients_and_refusals_do_not_mutate_physical_history() {
    for law in laws() {
        let held = law
            .bonded_point()
            .trial([-2e-4, 6e-4, 1e-4], 1e-8)
            .unwrap()
            .point;
        let zero_normal = law
            .bonded_point()
            .trial([0.0, 6e-4, 1e-4], 1e-8)
            .unwrap()
            .point;
        for (old, end) in [(held, held.jump_m()), (zero_normal, [0.0, 8e-4, 2e-4])] {
            let saved = old;
            let original = old
                .trial_over_time_with_damage_weighted_compression(end, 0.01, 1e-8)
                .unwrap();
            let result = old
                .trial_over_time_with_damage_weighted_compression_tangent(end, 0.01, 1e-8)
                .unwrap();
            assert_eq!((result.0, result.1, result.2), original);
            assert_eq!(result.4, None);
            assert_eq!(old, saved);
        }
        let unloaded = held.trial([-2e-4, 1e-4, 2e-5], 1e-8).unwrap().point;
        assert!(
            unloaded
                .trial_over_time_with_damage_weighted_compression_tangent(
                    unloaded.jump_m(),
                    0.01,
                    1e-8
                )
                .unwrap()
                .4
                .is_some()
        );

        let saved = held;
        for (end, dt, tolerance) in [
            ([0.03, 0.0, 0.0], 0.01, 1e-8),
            ([-2e-4, 8e-4, 1e-4], 0.0, 1e-8),
            ([-2e-4, f64::NAN, 1e-4], 0.01, 1e-8),
            ([-2e-4, 8e-4, 1e-4], 0.01, -1.0),
        ] {
            let original = held
                .trial_over_time_with_damage_weighted_compression(end, dt, tolerance)
                .unwrap_err();
            let result = held
                .trial_over_time_with_damage_weighted_compression_tangent(end, dt, tolerance)
                .unwrap_err();
            assert_eq!(result, original);
            assert_eq!(held, saved);
        }
    }
}

//! Local sliding/irreversible work in the original scalar fixed-load model.
//! Material values are synthetic; no empirical wear calibration is claimed.
use core::f64::consts::PI;
use fs_tribo::mindlin::surface::sliding::MindlinSurfaceAdvance;
use fs_tribo::mindlin::{MindlinDeresiewiczPoint, MindlinPointParameters};
use fs_tribo::{TangentialLaw, TangentialSolveConfig, TangentialStepInput};

fn model(mu: f64) -> MindlinDeresiewiczPoint {
    MindlinDeresiewiczPoint::new(MindlinPointParameters {
        normal_force_n: 10.0,
        contact_radius_m: 1e-3,
        shear_moduli_pa: [30e9, 20e9],
        poisson_ratios: [0.25, 0.3],
        friction_coefficient: mu,
        maximum_history_intervals: 64,
    })
    .unwrap()
}
fn input(du: f64) -> TangentialStepInput {
    TangentialStepInput {
        normal_n: 10.0,
        slip_increment_m: [du, 0.0],
        dt_s: 1e-3,
        contact_temperature_k: None,
        nominal_pressure_pa: None,
        compliance_m_per_n: 0.0,
        solve: TangentialSolveConfig {
            max_iterations: 64,
            maximum_force_n: 100.0,
            force_absolute_tolerance_n: 1e-13,
            force_relative_tolerance: 1e-13,
        },
    }
}
fn close(actual: f64, expected: f64, tol: f64) {
    assert!(
        (actual - expected).abs() <= tol,
        "actual={actual:e}, expected={expected:e}, tolerance={tol:e}"
    );
}
// Eq.4.58: macro travel minus elastic displacement in the virgin slip annulus.
fn virgin_sliding(radius: f64, fraction: f64) -> f64 {
    if fraction >= 1.0 {
        return fraction - 1.0 + 0.5 * radius * radius;
    }
    let excess = fraction - (1.0 - radius) * (1.0 + radius);
    if excess <= 0.0 {
        return 0.0;
    }
    let square = excess / (radius * radius);
    let t = square.sqrt();
    let shape = if t < 0.01 {
        t * t
            * t
            * (4.0 / 3.0 + square * (2.0 / 15.0 + square * (3.0 / 70.0 + square * 5.0 / 252.0)))
    } else {
        t * (1.0 - square).max(0.0).sqrt() - (1.0 - 2.0 * square) * t.asin()
    };
    radius * radius * shape / PI
}

#[test]
fn virgin_slip_is_zero_in_stick_and_matches_the_analytic_sliding_annulus() {
    let law = model(0.6);
    let a = law.parameters().contact_radius_m;
    let d = law.gross_slip_displacement_m();
    let before = law.initial_state();
    for fraction in [1e-8, 0.05, 0.5, 0.99, 1.0, 2.5] {
        for sign in [-1.0, 1.0] {
            let receipt = law
                .advance_with_surface_sliding(&before, &input(sign * fraction * d))
                .unwrap();
            let actual_fraction = receipt.tangential_step().slip_increment_m[0].abs() / d;
            for nominal in [0.0, 0.3, 0.7, 0.85, 0.99, 1.0] {
                let radius = nominal * a;
                let normalized = radius / a;
                let local = receipt.surface_sliding(radius).unwrap();
                let expected = d * virgin_sliding(normalized, actual_fraction);
                close(
                    local.sliding_distance_m,
                    expected,
                    2e-12 * expected + 1e-27 * d,
                );
                close(
                    local.signed_sliding_increment_m,
                    sign * expected,
                    2e-12 * expected + 1e-27 * d,
                );
                close(
                    local.dissipated_energy_j_per_m2,
                    0.6 * local.normal_pressure_pa * expected,
                    2e-12 * local.dissipated_energy_j_per_m2 + 1e-24,
                );
                if expected == 0.0 {
                    assert_eq!(local.sliding_distance_m, 0.0);
                }
                let old = law.surface_response(&before, radius).unwrap();
                let new = law
                    .surface_response(&receipt.tangential_step().state, radius)
                    .unwrap();
                close(
                    local.signed_sliding_increment_m,
                    receipt.tangential_step().slip_increment_m[0]
                        - (new.elastic_displacement_m - old.elastic_displacement_m),
                    4e-14 * d,
                );
            }
            if fraction < 1.0 {
                assert!(receipt.pressure_weighted_sliding_j() < 0.5 * 10.0 * fraction * d);
            }
        }
    }
}

#[test]
fn reversal_and_gross_sliding_have_their_actual_local_path_lengths() {
    let law = model(0.6);
    let a = law.parameters().contact_radius_m;
    let d = law.gross_slip_displacement_m();
    let virgin = law.initial_state();
    let peak = law.advance(&virgin, &input(0.8 * d)).unwrap().state;
    let reversed = law
        .advance_with_surface_sliding(&peak, &input(-0.6 * d))
        .unwrap();
    for normalized in [0.0, 0.3, 0.8, 0.9, 0.99, 1.0] {
        let local = reversed.surface_sliding(normalized * a).unwrap();
        let expected = 2.0 * d * virgin_sliding(normalized, 0.3);
        close(local.signed_sliding_increment_m, -expected, 3e-14 * d);
        close(local.sliding_distance_m, expected, 3e-14 * d);
    }
    // Once the entire patch slides in the existing direction, elastic motion
    // is unchanged: local sliding is the whole accepted increment everywhere.
    let saturated = law.advance(&virgin, &input(2.0 * d)).unwrap().state;
    let gross = law
        .advance_with_surface_sliding(&saturated, &input(0.4 * d))
        .unwrap();
    for normalized in [0.0, 0.2, 0.9, 1.0] {
        close(
            gross
                .surface_sliding(normalized * a)
                .unwrap()
                .sliding_distance_m,
            0.4 * d,
            2e-14 * d,
        );
    }
    close(
        gross.pressure_weighted_sliding_j(),
        10.0 * 0.4 * d,
        2e-14 * 10.0 * d,
    );
    // Full reversal first unloads and reloads the elastic field. It must not
    // charge the whole macro travel as plastic sliding.
    let full_reverse = law
        .advance_with_surface_sliding(&saturated, &input(-3.0 * d))
        .unwrap();
    for normalized in [0.0, 0.3, 0.9, 1.0] {
        close(
            full_reverse
                .surface_sliding(normalized * a)
                .unwrap()
                .sliding_distance_m,
            d * (1.0 + normalized * normalized),
            2e-14 * d,
        );
    }
    close(
        full_reverse.pressure_weighted_sliding_j(),
        1.4 * 10.0 * d,
        3e-14 * 10.0 * d,
    );
}

fn integrated(receipt: &MindlinSurfaceAdvance<'_>, a: f64, n: usize) -> [f64; 2] {
    let mut values = [0.0; 2];
    for i in 0..n {
        let radius = a * ((i as f64 + 0.5) / n as f64).sqrt();
        let local = receipt.surface_sliding(radius).unwrap();
        values[0] += local.dissipated_energy_j_per_m2;
        values[1] += local.normal_pressure_pa * local.sliding_distance_m;
    }
    values.map(|sum| PI * a * a * sum / n as f64)
}

#[test]
fn independent_area_integrals_converge_to_original_native_dissipation() {
    let law = model(0.6);
    let a = law.parameters().contact_radius_m;
    let d = law.gross_slip_displacement_m();
    let mut state = law.initial_state();
    for fraction in [0.7, -0.5, 0.2, 1.6, 0.3, -3.0] {
        let receipt = law
            .advance_with_surface_sliding(&state, &input(fraction * d))
            .unwrap();
        let native = receipt.tangential_step().dissipated_j;
        let errors = [512, 2048, 8192].map(|n| {
            let values = integrated(&receipt, a, n);
            close(values[0], 0.6 * values[1], 3e-14 * native);
            (values[0] - native).abs()
        });
        assert!(errors[2] < 0.2 * errors[0], "du/d={fraction}: {errors:?}");
        assert!(
            errors[2] < 3e-5 * native,
            "du/d={fraction}: {errors:?}, native={native:e}"
        );
        close(
            0.6 * receipt.pressure_weighted_sliding_j(),
            native,
            2e-15 * native,
        );
        state = receipt.into_tangential_step().state;
    }
    // Monotone subdivision preserves local path length and physical work.
    let direct = law
        .advance_with_surface_sliding(&state, &input(0.6 * d))
        .unwrap();
    let first = law
        .advance_with_surface_sliding(&state, &input(0.2 * d))
        .unwrap();
    let second = law
        .advance_with_surface_sliding(&first.tangential_step().state, &input(0.4 * d))
        .unwrap();
    for normalized in [0.0, 0.4, 0.8, 0.99, 1.0] {
        let radius = normalized * a;
        let all = direct.surface_sliding(radius).unwrap();
        let one = first.surface_sliding(radius).unwrap();
        let two = second.surface_sliding(radius).unwrap();
        close(
            all.sliding_distance_m,
            one.sliding_distance_m + two.sliding_distance_m,
            3e-14 * d,
        );
    }
}

#[test]
fn actual_compliant_motion_frictionless_limit_and_original_refusals_are_preserved() {
    let law = model(0.6);
    let before = law.initial_state();
    let d = law.gross_slip_displacement_m();
    let a = law.parameters().contact_radius_m;
    let mut request = input(0.8 * d);
    request.compliance_m_per_n = 2.0 / law.initial_stiffness_n_per_m();
    let native = law.advance(&before, &request).unwrap();
    let receipt = law.advance_with_surface_sliding(&before, &request).unwrap();
    assert_eq!(receipt.tangential_step(), &native);
    assert!(native.slip_increment_m[0] < request.slip_increment_m[0]);
    let local = receipt.surface_sliding(a).unwrap();
    close(
        local.sliding_distance_m,
        d * virgin_sliding(1.0, native.slip_increment_m[0] / d),
        3e-14 * d,
    );
    for radius in [-1.0, 1.01 * a, f64::NAN, f64::INFINITY] {
        assert!(receipt.surface_sliding(radius).is_err());
    }
    let free = model(0.0);
    assert!(
        free.advance_with_surface_sliding(&before, &input(d))
            .is_err()
    );
    request.normal_n = 11.0;
    assert!(law.advance_with_surface_sliding(&before, &request).is_err());
    let free_state = free.initial_state();
    let sliding = free
        .advance_with_surface_sliding(&free_state, &input(-d))
        .unwrap();
    assert_eq!(sliding.pressure_weighted_sliding_j(), 10.0 * d);
    for radius in [0.0, 0.5 * a, a] {
        let local = sliding.surface_sliding(radius).unwrap();
        assert_eq!(local.signed_sliding_increment_m, -d);
        assert_eq!(local.dissipated_energy_j_per_m2, 0.0);
    }
    let stationary = law
        .advance_with_surface_sliding(&before, &input(0.0))
        .unwrap();
    assert_eq!(stationary.pressure_weighted_sliding_j(), 0.0);
    assert_eq!(
        stationary
            .surface_sliding(0.5 * a)
            .unwrap()
            .sliding_distance_m,
        0.0
    );
    assert_eq!(before, law.initial_state());
}

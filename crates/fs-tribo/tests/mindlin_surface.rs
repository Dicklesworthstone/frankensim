//! Physical surface fields of the retained scalar Hertz--Mindlin continuum.
//! Synthetic material values verify the equations, not material calibration.
use core::f64::consts::PI;
use fs_tribo::mindlin::loading::{HertzMindlinContact, HertzMindlinParameters};
use fs_tribo::mindlin::{MindlinDeresiewiczPoint, MindlinPointParameters};
use fs_tribo::{TangentialLaw, TangentialSolveConfig, TangentialStepInput};

fn model() -> MindlinDeresiewiczPoint {
    MindlinDeresiewiczPoint::new(MindlinPointParameters {
        normal_force_n: 10.0,
        contact_radius_m: 1e-3,
        shear_moduli_pa: [30e9, 20e9],
        poisson_ratios: [0.25, 0.3],
        friction_coefficient: 0.6,
        maximum_history_intervals: 64,
    })
    .unwrap()
}
fn input(normal_n: f64, du: f64) -> TangentialStepInput {
    TangentialStepInput {
        normal_n,
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
fn hertz_profile(radius: f64, patch: f64) -> f64 {
    ((patch - radius) * (patch + radius)).max(0.0).sqrt()
}
// Independent closed form of a parabolic MDR field, Eq. 4.58. The virgin
// partial-slip displacement is the difference of two such disk fields.
fn disk_displacement(radius: f64, patch: f64) -> f64 {
    if radius <= patch {
        patch * patch - 0.5 * radius * radius
    } else {
        ((2.0 * patch * patch - radius * radius) * (patch / radius).asin()
            + patch * (radius * radius - patch * patch).sqrt())
            / PI
    }
}
fn virgin_displacement(radius: f64, fraction: f64) -> f64 {
    let fraction = fraction.min(1.0);
    let c = (1.0 - fraction).sqrt();
    if radius <= c {
        fraction
    } else if radius <= 1.0 {
        // Eq. 4.58 with asin(c/r)=pi/2-asin(sqrt(1-c^2/r^2)). This
        // avoids subtracting nearly equal disk fields at very small slip.
        let square = (fraction - (1.0 - radius) * (1.0 + radius)) / (radius * radius);
        let t = square.max(0.0).sqrt();
        fraction
            + radius * radius / PI
                * ((1.0 - 2.0 * square) * t.asin() - t * (1.0 - square).max(0.0).sqrt())
    } else {
        disk_displacement(radius, 1.0) - disk_displacement(radius, c)
    }
}

#[test]
fn virgin_fields_match_cattaneo_pressure_superposition_and_surface_motion() {
    let law = model();
    let p = law.parameters();
    let d = law.gross_slip_displacement_m();
    let p0 = 1.5 * p.normal_force_n / (PI * p.contact_radius_m.powi(2));
    for fraction in [1e-8_f64, 0.02, 0.4, 0.999, 1.0, 3.0] {
        for sign in [-1.0, 1.0] {
            let increment = sign * fraction * d;
            let fraction = increment.abs() / d;
            let c = (1.0 - fraction.min(1.0)).sqrt();
            let state = law
                .advance(&law.initial_state(), &input(p.normal_force_n, increment))
                .unwrap()
                .state;
            for radius in [0.0, 1e-150, 0.25, 0.7, c, 0.999999999, 1.0, 1.01, 3.0] {
                let query_radius = radius * p.contact_radius_m;
                let radius = query_radius / p.contact_radius_m;
                let field = law.surface_response(&state, query_radius).unwrap();
                let pressure = p0 * hertz_profile(radius, 1.0);
                // Use c^2=1-u/d directly: rounding sqrt(c^2) and squaring it
                // again moves the singular stick boundary by a floating ULP.
                let inner = ((1.0 - radius) * (1.0 + radius) - fraction.min(1.0))
                    .max(0.0)
                    .sqrt();
                let traction =
                    sign * p.friction_coefficient * p0 * (hertz_profile(radius, 1.0) - inner);
                let displacement = sign * d * virgin_displacement(radius, fraction);
                close(field.normal_pressure_pa, pressure, 2e-14 * p0);
                close(field.tangential_traction_pa, traction, 3e-14 * p0);
                assert!(
                    (field.elastic_displacement_m - displacement).abs() <= 3e-14 * d,
                    "fraction={fraction}, radius={radius}, sign={sign}: actual={:e}, expected={displacement:e}",
                    field.elastic_displacement_m
                );
                if fraction <= 1.0 && radius <= c {
                    close(field.elastic_displacement_m, sign * fraction * d, 3e-14 * d);
                }
            }
            let radius = 1e12 * p.contact_radius_m;
            let far = law.surface_response(&state, radius).unwrap();
            let g_mdr = law.initial_stiffness_n_per_m() / (2.0 * p.contact_radius_m);
            let asymptote = state.endpoint_force_n() / (PI * g_mdr * radius);
            close(
                far.elastic_displacement_m,
                asymptote,
                1e-12 * asymptote.abs(),
            );
            assert_eq!(far.normal_pressure_pa, 0.0);
            assert_eq!(far.tangential_traction_pa, 0.0);
        }
    }
}

#[test]
fn reversal_retains_opposite_signed_outer_traction_and_exact_superposition() {
    let law = model();
    let p = law.parameters();
    let d = law.gross_slip_displacement_m();
    let virgin = law.initial_state();
    let peak = law
        .advance(&virgin, &input(p.normal_force_n, 0.8 * d))
        .unwrap()
        .state;
    let reversed = law
        .advance(&peak, &input(p.normal_force_n, -0.6 * d))
        .unwrap()
        .state;
    let correction = law
        .advance(&virgin, &input(p.normal_force_n, 0.3 * d))
        .unwrap()
        .state;
    let p0 = 1.5 * p.normal_force_n / (PI * p.contact_radius_m.powi(2));
    for index in 0..=200 {
        let radius = 1.5 * p.contact_radius_m * index as f64 / 200.0;
        let actual = law.surface_response(&reversed, radius).unwrap();
        let a = law.surface_response(&peak, radius).unwrap();
        let b = law.surface_response(&correction, radius).unwrap();
        close(
            actual.tangential_traction_pa,
            a.tangential_traction_pa - 2.0 * b.tangential_traction_pa,
            2e-14 * p0,
        );
        close(
            actual.elastic_displacement_m,
            a.elastic_displacement_m - 2.0 * b.elastic_displacement_m,
            2e-14 * d,
        );
        assert!(
            actual.tangential_traction_pa.abs()
                <= p.friction_coefficient * actual.normal_pressure_pa + 2e-14 * p0
        );
    }
    assert!(
        law.surface_response(&reversed, 0.0)
            .unwrap()
            .tangential_traction_pa
            > 0.0
    );
    assert!(
        law.surface_response(&reversed, 0.99 * p.contact_radius_m)
            .unwrap()
            .tangential_traction_pa
            < 0.0
    );
    let saved = reversed.clone();
    let errors: Vec<_> = [256, 1024, 4096]
        .map(|n| {
            // Independent equal-area midpoint integration of the returned field.
            let mut force = 0.0;
            for i in 0..n {
                let radius = p.contact_radius_m * ((i as f64 + 0.5) / n as f64).sqrt();
                force += law
                    .surface_response(&reversed, radius)
                    .unwrap()
                    .tangential_traction_pa;
            }
            (PI * p.contact_radius_m.powi(2) * force / n as f64 - reversed.endpoint_force_n()).abs()
        })
        .into_iter()
        .collect();
    assert!(errors[2] < 0.2 * errors[0], "{errors:?}");
    assert!(errors[2] < 5e-5 * law.gross_slip_capacity_n(), "{errors:?}");
    assert_eq!(saved, reversed);
}

#[test]
fn simultaneous_normal_shear_growth_and_unloading_reconstruct_current_geometry() {
    let law = HertzMindlinContact::new(HertzMindlinParameters {
        effective_radius_m: 0.02,
        shear_moduli_pa: [30e9, 20e9],
        poisson_ratios: [0.25, 0.3],
        friction_coefficient: 0.6,
        maximum_history_intervals: 64,
    })
    .unwrap();
    let open = law.initial_state();
    let normal = law.normal_step(&open, 10.0).unwrap().state;
    let fixed = MindlinDeresiewiczPoint::new(MindlinPointParameters {
        normal_force_n: 10.0,
        contact_radius_m: normal.contact_radius_m(),
        shear_moduli_pa: law.parameters().shear_moduli_pa,
        poisson_ratios: law.parameters().poisson_ratios,
        friction_coefficient: law.parameters().friction_coefficient,
        maximum_history_intervals: 64,
    })
    .unwrap();
    let d = fixed.gross_slip_displacement_m();
    let grown = law
        .advance_prescribed_path(&open, &input(10.0, 0.4 * d), &mut || false)
        .unwrap()
        .tangential
        .state;
    let a = grown.contact_radius_m();
    for normalized in [0.0, 0.2, 0.9, 1.0, 2.0] {
        let f = law.surface_response(&grown, a * normalized).unwrap();
        close(
            f.tangential_traction_pa,
            0.4 * 0.6 * f.normal_pressure_pa,
            2e-14 * 10.0 / (a * a),
        );
        close(
            f.elastic_displacement_m,
            0.4 * d * disk_displacement(normalized, 1.0),
            3e-14 * d,
        );
    }
    let mut state = grown;
    for (normal_n, du) in [(3.0, -0.15 * d), (14.0, 0.12 * d), (7.0, -0.4 * d)] {
        state = law
            .advance_prescribed_path(&state, &input(normal_n, du), &mut || false)
            .unwrap()
            .tangential
            .state;
        let a = state.contact_radius_m();
        for i in 0..=64 {
            let rho = a * i as f64 / 64.0;
            let f = law.surface_response(&state, rho).unwrap();
            assert!(
                f.tangential_traction_pa.abs()
                    <= 0.6 * f.normal_pressure_pa + 1e-13 * normal_n / (a * a)
            );
            assert!(f.elastic_displacement_m.is_finite());
        }
    }
    let separated = law.normal_step(&state, 0.0).unwrap().state;
    let recontact = law.normal_step(&separated, 10.0).unwrap().state;
    for radius in [0.0, a, 2.0 * a] {
        assert_eq!(
            law.surface_response(&separated, radius).unwrap(),
            Default::default()
        );
        let f = law.surface_response(&recontact, radius).unwrap();
        assert_eq!(f.tangential_traction_pa, 0.0);
        assert_eq!(f.elastic_displacement_m, 0.0);
    }
}

#[test]
fn original_identity_and_domain_are_enforced_without_mutating_history() {
    let law = model();
    let state = law.initial_state();
    let mut parameters = law.parameters();
    parameters.shear_moduli_pa[0] *= 2.0;
    let foreign = MindlinDeresiewiczPoint::new(parameters).unwrap();
    assert!(foreign.surface_response(&state, 0.0).is_err());
    for radius in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(law.surface_response(&state, radius).is_err());
    }
    parameters = law.parameters();
    parameters.friction_coefficient = 0.0;
    let free = MindlinDeresiewiczPoint::new(parameters).unwrap();
    let moved = free
        .advance(
            &free.initial_state(),
            &input(parameters.normal_force_n, 1.0),
        )
        .unwrap()
        .state;
    for radius in [
        0.0,
        0.5 * parameters.contact_radius_m,
        parameters.contact_radius_m,
        f64::MAX,
    ] {
        let f = free.surface_response(&moved, radius).unwrap();
        assert_eq!(f.tangential_traction_pa, 0.0);
        assert_eq!(f.elastic_displacement_m, 0.0);
    }
    assert_eq!(state, law.initial_state());
}

//! G0/G1/G3: an explicit steady sliding law, not thermal rosin validation.
use fs_tribo::{FrictionLaw, FrictionRegime, TangentialLaw, TangentialSolveConfig, TangentialStepInput};

fn law() -> FrictionLaw {
    FrictionLaw::Hyperbolic { static_mu: 0.8, kinetic_mu: 0.3, characteristic_speed: 0.04 }
}
fn input(distance: f64, compliance: f64, normal: f64, dt: f64) -> TangentialStepInput {
    TangentialStepInput {
        normal_n: normal, slip_increment_m: [distance, 0.0], dt_s: dt,
        contact_temperature_k: None, nominal_pressure_pa: None, compliance_m_per_n: compliance,
        solve: TangentialSolveConfig { max_iterations: 96, maximum_force_n: 100.0,
            force_absolute_tolerance_n: 1e-12, force_relative_tolerance: 1e-11 },
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn g1_hyperbolic_tail_and_independent_numerical_regularization() {
    let law = law();
    assert_eq!(law.static_coefficient().unwrap(), 0.8);
    assert_eq!(law.zero_speed_kinetic_coefficient().unwrap(), 0.8);
    for speed in [1e-9, 0.004, 0.04, 0.4, 4.0, 400.0] {
        let expected = 0.3 + 0.02 / (speed + 0.04);
        close(law.kinetic_coefficient(speed).unwrap(), expected, 2e-15);
        for regularization in [1e-12, 1e-10] {
            let force = law.regularized_traction_1d(speed, 2.0, regularization).unwrap();
            close(force, -2.0*expected, 4e-15);
            assert_eq!(force, -law.regularized_traction_1d(-speed, 2.0, regularization).unwrap());
            assert!(-force*speed >= 0.0);
        }
    }
    let large = FrictionLaw::Hyperbolic { static_mu: 0.8, kinetic_mu: 0.3,
        characteristic_speed: f64::MAX };
    close(large.kinetic_coefficient(f64::MAX).unwrap(), 0.55, 2e-15);
    let tiny = FrictionLaw::Hyperbolic { static_mu: 0.8, kinetic_mu: 0.3,
        characteristic_speed: f64::from_bits(1) };
    assert_eq!(tiny.kinetic_coefficient(f64::MAX).unwrap(), 0.3);
}

#[test]
fn g1_implicit_sliding_matches_independent_quadratic_and_actual_work() {
    let law = law();
    for dt in [1.0/48000.0, 0.001] {
        for compliance in [1e-8, 1e-6, 1e-4] {
            for normal in [0.1, 1.0, 2.0] {
                let distance = 1.3*compliance*normal*0.8 + dt*0.1;
                // Eliminate velocity from F=N*(mu_d+K/(v+v0)) and
                // dt*v=distance-compliance*F. Select the smaller quadratic root.
                let b: f64 = distance + dt*0.04 + compliance*normal*0.3;
                let c = normal*(0.3*distance + 0.8*dt*0.04);
                let reference = 2.0*c/(b + (b*b-4.0*compliance*c).sqrt());
                for sign in [-1.0, 1.0] {
                    let initial = law.initial_state();
                    let step = law.advance(&initial, &input(sign*distance,compliance,normal,dt)).unwrap();
                    assert_eq!(initial, law.initial_state(), "advance returns a candidate");
                    assert_eq!(step.regime, FrictionRegime::Sliding);
                    assert_eq!(step.iterations, 0, "quadratic sliding needs no bisection");
                    close(step.traction_n[0], sign*reference, 5e-11);
                    close(step.slip_increment_m[0], sign*distance-compliance*step.traction_n[0], 1e-20);
                    assert_eq!(step.dissipated_j, step.traction_n[0]*step.slip_increment_m[0]);
                    assert!(step.dissipated_j >= 0.0);
                    assert_eq!(step.closure_residual_j, 0.0);
                }
            }
        }
    }
}

#[test]
fn g1_direct_hyperbolic_root_preserves_the_static_selection_and_force_budget() {
    let law = law();
    // Both equality and strict static feasibility must remain with the owner's
    // static test, even when the weakening load line has another sliding root.
    for distance in [0.1, 0.5, 0.8] {
        assert_eq!(
            law.direct_sliding_traction_n(1.0, distance, 1.0, 0.01).unwrap(),
            None,
        );
        let step = law.advance(&law.initial_state(), &input(distance, 1.0, 1.0, 0.01)).unwrap();
        assert_eq!(step.regime, FrictionRegime::Sticking);
        assert_eq!(step.dissipated_j, 0.0);
    }
    let mut x = input(0.001, 0.001, 1.0, 0.001);
    x.solve.max_iterations = 1;
    let initial = law.initial_state();
    let step = law.advance(&initial, &x).unwrap();
    assert_eq!(step.regime, FrictionRegime::Sliding);
    assert_eq!(step.iterations, 0);
    x.solve.maximum_force_n = 0.2;
    assert!(law.advance(&initial, &x).is_err());
    assert_eq!(initial, law.initial_state());
}

#[test]
fn g3_direct_hyperbolic_root_survives_dimensional_scaling_and_near_stick() {
    let law = law();
    // The raw quadratic squares lengths and multiplies length by load. Its
    // coefficients overflow/underflow on these otherwise representable cases.
    for distance_scale in [1e-200, 1.0, 1e200] {
        for force_scale in [1e-100, 1.0, 1e100] {
            let distance = 0.003 * distance_scale;
            let compliance = 0.001 * distance_scale / force_scale;
            let dt = 0.001 * distance_scale;
            let reference = law.direct_sliding_traction_n(1.0, 0.003, 0.001, 0.001)
                .unwrap().unwrap();
            let force = law.direct_sliding_traction_n(force_scale, distance, compliance, dt)
                .unwrap().unwrap();
            close(force / force_scale, reference, 3e-15);
            let slip = distance - compliance * force;
            assert!(slip > 0.0);
            let expected = force_scale * law.kinetic_coefficient(slip / dt).unwrap();
            close(force / force_scale, expected / force_scale, 3e-15);
        }
    }
    // At the limiting static state the slip quadratic has a double root when
    // d/dt=0.064 m/s. Nearby roots must satisfy the original law without a
    // negative discriminant from subtracting nearly equal squared terms.
    for excess in [1e-12, 1e-9, 1e-6] {
        for released_speed in [0.0064, 0.064, 0.64] {
            let distance = 0.8 * (1.0 + excess);
            let dt = distance / released_speed;
            let force = law.direct_sliding_traction_n(1.0, distance, 1.0, dt)
                .unwrap().unwrap();
            assert!(force < distance);
            let expected = law.kinetic_coefficient((distance - force) / dt).unwrap();
            close(force, expected, 5e-15);
        }
    }
}

#[test]
fn g1_direct_hyperbolic_root_matches_an_independent_bisection() {
    let law = law();
    for normal in [0.01, 0.7, 3.0] {
        for compliance in [1e-8, 1e-4, 0.1] {
            for speed in [1e-5, 0.04, 10.0] {
                let dt = 1.0 / 48_000.0;
                let distance = compliance * normal * 0.8 * 1.01 + dt * speed;
                let mut lo = 0.0;
                let mut hi = normal * 0.8;
                for _ in 0..100 {
                    let mid = f64::midpoint(lo, hi);
                    let actual_speed = (distance - compliance * mid) / dt;
                    let expected = normal * (0.3 + 0.02 / (actual_speed + 0.04));
                    if mid < expected { lo = mid; } else { hi = mid; }
                }
                let force = law.direct_sliding_traction_n(normal, distance, compliance, dt)
                    .unwrap().unwrap();
                close(force, f64::midpoint(lo, hi), 5e-15 * normal);
            }
        }
    }
}

#[test]
fn g0_sticking_zero_loss_and_refusal_preserve_constitutive_semantics() {
    let law = law();
    let state = law.initial_state();
    let step = law.advance(&state, &input(1e-8,1e-6,1.0,0.001)).unwrap();
    assert_eq!(step.regime, FrictionRegime::Sticking);
    close(step.traction_n[0], 0.01, 1e-15);
    assert_eq!(step.slip_increment_m, [0.0,0.0]);
    assert_eq!(step.dissipated_j, 0.0);
    let zero = FrictionLaw::Hyperbolic { static_mu: 0.0, kinetic_mu: 0.0, characteristic_speed: 0.04 };
    assert_eq!(zero.advance(&zero.initial_state(), &input(1e-5,1e-6,1.0,0.001)).unwrap().traction_n, [0.0,0.0]);
    for (static_mu, kinetic_mu, characteristic_speed) in [
        (-1.0,0.3,0.04), (0.8,0.9,0.04), (0.8,-0.1,0.04),
        (0.8,0.3,0.0), (0.8,0.3,f64::NAN), (f64::INFINITY,0.3,0.04),
    ] {
        let bad = FrictionLaw::Hyperbolic { static_mu, kinetic_mu, characteristic_speed };
        assert!(bad.validate().is_err());
        assert!(bad.static_coefficient().is_err());
        assert!(bad.zero_speed_kinetic_coefficient().is_err());
    }
    let mut bounded = input(1e-3,1e-6,1.0,0.001);
    bounded.solve.maximum_force_n = 0.01;
    assert!(law.advance(&state,&bounded).is_err());
    assert_eq!(state,law.initial_state());
}

//! Synthetic constitutive checks of real-area adhesive held-effort friction.
use fs_tribo::adhesive::{AdhesiveContactState, AdhesiveFrictionLaw, AdhesiveTangentialContact};
use fs_tribo::{
    FrictionLaw, FrictionRegime, TangentialLaw, TangentialSolveConfig, TangentialStepInput,
};
fn input(normal_n: f64, slip: [f64; 2], compliance: f64) -> TangentialStepInput {
    TangentialStepInput {
        normal_n,
        slip_increment_m: slip,
        dt_s: 1.0,
        contact_temperature_k: None,
        nominal_pressure_pa: None,
        compliance_m_per_n: compliance,
        solve: TangentialSolveConfig {
            max_iterations: 128,
            maximum_force_n: 100.0,
            force_absolute_tolerance_n: 1e-13,
            force_relative_tolerance: 1e-13,
        },
    }
}
fn contact(
    law: FrictionLaw,
    area: f64,
    normal: f64,
    ts: f64,
    tk: f64,
) -> AdhesiveTangentialContact {
    AdhesiveTangentialContact::new(
        AdhesiveFrictionLaw::new(law, ts, tk).unwrap(),
        AdhesiveContactState::new(area, normal).unwrap(),
    )
    .unwrap()
}
fn zero_pressure() -> FrictionLaw {
    FrictionLaw::Coulomb {
        static_mu: 0.0,
        kinetic_mu: 0.0,
    }
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 5e-12 * a.abs().max(b.abs()).max(1e-20),
        "{a} != {b}"
    );
}
#[test]
fn zero_repulsive_load_sticks_and_slides_with_vector_work() {
    let law = contact(zero_pressure(), 2e-6, 0.0, 3e6, 2e6);
    let initial = law.initial_state();
    let stuck = law.advance(&initial, &input(0.0, [0.3, 0.4], 0.1)).unwrap();
    assert_eq!(stuck.regime, FrictionRegime::Sticking);
    close(stuck.traction_n[0], 3.0);
    close(stuck.traction_n[1], 4.0);
    assert_eq!(stuck.dissipated_j, 0.0);
    assert_eq!(stuck.slip_increment_m, [0.0; 2]);
    let slide = law
        .advance(&stuck.state, &input(0.0, [0.6, 0.8], 0.1))
        .unwrap();
    assert_eq!(slide.regime, FrictionRegime::Sliding);
    close(slide.traction_n[0], 2.4);
    close(slide.traction_n[1], 3.2);
    close(slide.slip_increment_m[0], 0.36);
    close(slide.slip_increment_m[1], 0.48);
    close(slide.dissipated_j, 2.4);
    assert_eq!(slide.stored_change_j, 0.0);
    assert_eq!(slide.closure_residual_j, 0.0);
    assert_ne!(
        slide.dissipated_j,
        slide.traction_n[0] * 0.6 + slide.traction_n[1] * 0.8
    );
    let reverse = law
        .advance(&initial, &input(0.0, [-0.6, -0.8], 0.1))
        .unwrap();
    close(reverse.dissipated_j, slide.dissipated_j);
    close(reverse.traction_n[0], -slide.traction_n[0]);
}
#[test]
fn affine_pressure_and_junction_shear_use_the_physical_load_line() {
    for slope in [0.4, -0.6] {
        let law = contact(
            FrictionLaw::VelocityDependent {
                static_mu: 0.4,
                mu_zero: 0.3,
                slope_per_speed: slope,
            },
            1.0,
            1.0,
            0.5,
            0.5,
        );
        let x = input(1.0, [0.6, 0.0], 0.3);
        let out = law.advance(&law.initial_state(), &x).unwrap();
        assert_eq!(out.regime, FrictionRegime::Sliding);
        assert_eq!(out.iterations, 0);
        let expected = (0.5 + 0.3 + slope * 0.6) / (1.0 + slope * 0.3);
        close(out.traction_n[0], expected);
        close(
            out.traction_n[0],
            0.5 + (0.3 + slope * out.slip_increment_m[0]).max(0.0),
        );
    }
    let law = contact(
        FrictionLaw::VelocityDependent {
            static_mu: 0.4,
            mu_zero: 0.3,
            slope_per_speed: -0.6,
        },
        1.0,
        1.0,
        0.5,
        0.5,
    );
    let out = law
        .advance(&law.initial_state(), &input(1.0, [2.0, 0.0], 0.3))
        .unwrap();
    close(out.traction_n[0], 0.5);
    assert_eq!(out.iterations, 0);
}
#[test]
fn zero_junction_strength_recovers_every_original_pressure_law() {
    let laws = [
        FrictionLaw::Coulomb {
            static_mu: 0.8,
            kinetic_mu: 0.3,
        },
        FrictionLaw::VelocityDependent {
            static_mu: 0.8,
            mu_zero: 0.4,
            slope_per_speed: 0.1,
        },
        FrictionLaw::Stribeck {
            static_mu: 0.8,
            kinetic_mu: 0.3,
            characteristic_speed: 0.05,
            viscous_per_speed: 0.01,
        },
        FrictionLaw::Hyperbolic {
            static_mu: 0.8,
            kinetic_mu: 0.3,
            characteristic_speed: 0.05,
        },
    ];
    for dry in laws {
        let law = contact(dry.clone(), 1.0, 2.0, 0.0, 0.0);
        for du in [[0.0, 0.0], [0.001, -0.002], [0.3, 0.4], [-0.8, 0.0]] {
            for c in [0.0, 0.01, 0.8] {
                let x = input(2.0, du, c);
                assert_eq!(
                    law.advance(&law.initial_state(), &x).unwrap(),
                    dry.advance(&dry.initial_state(), &x).unwrap()
                );
            }
        }
    }
}
#[test]
fn detached_patch_has_no_traction_or_heat_and_failed_steps_keep_history() {
    let attached = contact(zero_pressure(), 2e-6, 0.0, 3e6, 2e6);
    let prior = attached
        .advance(&attached.initial_state(), &input(0.0, [0.6, 0.8], 0.1))
        .unwrap()
        .state;
    let before = prior.clone();
    let detached = contact(zero_pressure(), 0.0, 0.0, 3e6, 2e6);
    let out = detached
        .advance(&prior, &input(0.0, [0.2, 0.0], 0.1))
        .unwrap();
    assert_eq!(out.traction_n, [0.0; 2]);
    assert_eq!(out.dissipated_j, 0.0);
    assert_eq!(out.slip_increment_m, [0.2, 0.0]);
    let mut wrong = input(1.0, [0.2, 0.0], 0.1);
    assert!(attached.advance(&prior, &wrong).is_err());
    wrong.normal_n = 0.0;
    wrong.dt_s = 0.0;
    assert!(attached.advance(&prior, &wrong).is_err());
    wrong.dt_s = 1.0;
    wrong.solve.maximum_force_n = 0.1;
    assert!(attached.advance(&prior, &wrong).is_err());
    assert_eq!(prior, before);
}
#[test]
fn nonlinear_pressure_capacity_closes_at_actual_compliant_slip() {
    let law = contact(
        FrictionLaw::Stribeck {
            static_mu: 0.8,
            kinetic_mu: 0.3,
            characteristic_speed: 0.1,
            viscous_per_speed: 0.1,
        },
        0.5,
        2.0,
        3.0,
        2.0,
    );
    let x = input(2.0, [0.6, 0.8], 0.1);
    let out = law.advance(&law.initial_state(), &x).unwrap();
    let slip = out.slip_increment_m[0].hypot(out.slip_increment_m[1]);
    let f = out.traction_n[0].hypot(out.traction_n[1]);
    close(
        f,
        law.law()
            .sliding_capacity_n(law.contact(), slip / x.dt_s)
            .unwrap(),
    );
    close(out.dissipated_j, f * slip);
}

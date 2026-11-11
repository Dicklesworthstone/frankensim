//! Original nonlinear wall law with a prescribed length, not constant beta.

use super::*;

#[test]
fn original_carreau_stress_selects_the_wall_jump_and_separate_loss() {
    let (b, wall, dt) = (0.4, -0.15, 0.05);
    let drive = FilmDrive {
        acceleration_m_s2: 0.8,
        wall_velocity_m_s: wall,
        surface_shear_pa: 0.15,
    };
    let tau_wall = drive.surface_shear_pa + SPEC.density_kg_m3 * drive.acceleration_m_s2;
    // Parent oracle uses an independent explicit Carreau formula and bisection.
    let wall_rate = reference_rate(tau_wall);
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let dy = 1.0 / n as f64;
        let mut u = wall + (b + 0.5 * dy) * wall_rate;
        let mut initial = vec![u];
        for face in 1..n {
            u += dy
                * reference_rate(
                    drive.surface_shear_pa + drive.acceleration_m_s2 * (1.0 - face as f64 * dy),
                );
            initial.push(u);
        }
        let mut film =
            TransientCarreauFilm::from_profile(law(0.5), 300.0, SPEC, initial.clone(), 0.0)
                .unwrap();
        film.set_wall_slip_length_m(b).unwrap();
        let mut heat = vec![f64::NAN; n];
        let report = film
            .step_with_heat(
                dt,
                drive,
                LIMITS,
                TransientCarreauControls::default(),
                &mut heat,
                || false,
            )
            .unwrap();
        for (&u, &v) in film.velocities_m_s().iter().zip(&initial) {
            near(u, v, 2e-12);
        }
        near(
            report.flow.midpoint_wall_fluid_velocity_m_s,
            wall + b * wall_rate,
            2e-13,
        );
        near(report.flow.midpoint_wall_shear_pa, tau_wall, 2e-12);
        near(
            report.flow.wall_slip_heat_j_m2,
            dt * tau_wall * b * wall_rate,
            2e-13,
        );
        near(heat.iter().sum(), report.flow.viscous_heat_j_m2, 2e-13);
        near(report.flow.work_residual_j_m2, 0.0, 2e-12);
        // A constant beta based on mu(0) would be a different boundary law: its slip speed is b tau_wall / mu(0) with
        // mu(0) = 2 Pa s, where the nonlinear law's is b wall_rate. With these inputs (printed by a run, n = 8):
        // tau_wall = 0.15 + 1 x 0.8 = 0.95 Pa, wall_rate = 0.5680108 1/s, so the nonlinear slip is 0.4 x 0.5680108 = 0.22720
        // m/s against 0.4 x 0.95 / 2 = 0.19 m/s, a difference of 0.0372 m/s, 16.4 % of the nonlinear slip. The threshold is
        // derived from that law, not a literal: the two boundary laws must differ by at least a tenth of the nonlinear slip.
        let nonlinear_slip = b * wall_rate;
        let constant_beta_slip = b * tau_wall / 2.0;
        assert!(
            (nonlinear_slip - constant_beta_slip).abs() > 0.1 * nonlinear_slip,
            "{nonlinear_slip} against {constant_beta_slip}"
        );
        let error = film
            .velocities_m_s()
            .iter()
            .enumerate()
            .map(|(i, &u)| {
                (u - wall - b * wall_rate - reference_velocity((i as f64 + 0.5) * dy)).abs()
            })
            .fold(0.0_f64, f64::max);
        errors.push(error);
    }
    assert!(errors[0] / errors[1] > 3.7, "{errors:?}");
    assert!(errors[1] / errors[2] > 3.7, "{errors:?}");
    assert!(errors[2] < 2e-4, "{errors:?}");
}

#[test]
fn slip_uses_original_thermal_preparation_and_newtonian_limit_without_resetting_history() {
    for b in [0.0, 0.2, 1.0] {
        let source = law(1.0);
        let mut nonlinear = TransientCarreauFilm::new(source, 300.0, SPEC, 16).unwrap();
        nonlinear.set_wall_slip_length_m(b).unwrap();
        nonlinear.set_temperature(330.0, || false).unwrap();
        assert_eq!(nonlinear.wall_slip_length_m(), b);
        assert_eq!(nonlinear.material().law(), &source);
        let mu = source
            .evaluate(330.0, 0.0)
            .unwrap()
            .dynamic_viscosity_pa_s();
        let mut linear = TransientFilm::new(
            FilmMaterial {
                density_kg_m3: 1.0,
                kinematic_viscosity_m2_s: mu,
                surface_tension_n_m: 0.0,
            },
            1.0,
            16,
        )
        .unwrap();
        linear.set_wall_slip_length_m(b).unwrap();
        for acceleration in [0.4, -0.7, 0.0] {
            let drive = FilmDrive {
                acceleration_m_s2: acceleration,
                wall_velocity_m_s: 0.03,
                surface_shear_pa: -0.1,
            };
            let a = nonlinear
                .step(
                    0.02,
                    drive,
                    LIMITS,
                    TransientCarreauControls::default(),
                    || false,
                )
                .unwrap();
            let z = linear.step(0.02, drive, LIMITS, || false).unwrap();
            for (&u, &v) in nonlinear
                .velocities_m_s()
                .iter()
                .zip(linear.velocities_m_s())
            {
                near(u, v, 3e-12);
            }
            near(a.flow.viscous_heat_j_m2, z.viscous_heat_j_m2, 1e-12);
            near(a.flow.wall_slip_heat_j_m2, z.wall_slip_heat_j_m2, 1e-12);
            near(
                a.flow.midpoint_wall_fluid_velocity_m_s,
                z.midpoint_wall_fluid_velocity_m_s,
                1e-12,
            );
            near(a.flow.work_residual_j_m2, 0.0, 2e-12);
        }
    }
}

#[test]
fn zero_slip_parity_and_nonzero_slip_refusal_retain_source_state_and_heat() {
    let source = law(0.45);
    let make = || TransientCarreauFilm::new(source, 300.0, SPEC, 16).unwrap();
    let mut a = make();
    let mut b = make();
    b.set_wall_slip_length_m(0.0).unwrap();
    let drive = FilmDrive {
        acceleration_m_s2: 0.8,
        wall_velocity_m_s: -0.1,
        surface_shear_pa: -0.2,
    };
    let controls = TransientCarreauControls::default();
    let first = a.step(0.04, drive, LIMITS, controls, || false).unwrap();
    assert_eq!(
        first,
        b.step(0.04, drive, LIMITS, controls, || false).unwrap()
    );
    for (&u, &v) in a.velocities_m_s().iter().zip(b.velocities_m_s()) {
        assert_eq!(u.to_bits(), v.to_bits());
    }
    for owner in [&mut a, &mut b] {
        owner.set_wall_slip_length_m(0.3).unwrap();
    }
    let mut expected_heat = vec![0.0; 16];
    let mut polls = 0;
    let expected = a
        .step_with_heat(0.04, drive, LIMITS, controls, &mut expected_heat, || {
            polls += 1;
            false
        })
        .unwrap();
    let saved = b.velocities_m_s().to_vec();
    let saved_material = *b.material();
    let time = b.accepted_time_s();
    let mut heat = vec![-5.0; 16];
    assert!(
        b.step_with_heat(
            0.04,
            drive,
            LIMITS,
            TransientCarreauControls {
                maximum_material_evaluations: 1,
                ..controls
            },
            &mut heat,
            || false
        )
        .is_err()
    );
    let mut seen = 0;
    assert!(matches!(
        b.step_with_heat(0.04, drive, LIMITS, controls, &mut heat, || {
            seen += 1;
            seen == polls
        }),
        Err(TransientCarreauError::Flow(TransientFilmError::Cancelled))
    ));
    assert_eq!(b.velocities_m_s(), saved);
    assert_eq!(b.accepted_time_s(), time);
    assert_eq!(*b.material(), saved_material);
    assert_eq!(b.wall_slip_length_m(), 0.3);
    assert!(heat.iter().all(|&v| v == -5.0));
    assert_eq!(
        b.step_with_heat(0.04, drive, LIMITS, controls, &mut heat, || false)
            .unwrap(),
        expected
    );
    assert_eq!(heat, expected_heat);
}

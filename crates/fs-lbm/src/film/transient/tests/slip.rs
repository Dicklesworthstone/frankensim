//! Original momentum with finite Navier length: independent steady/2-cell
//! references, signed source work, bulk/interface heat and transaction laws.

use super::*;

#[test]
fn zero_length_keeps_native_bits_and_invalid_rebinding_is_atomic() {
    let mut original = TransientFilm::new(MATERIAL, 1.0, 17).unwrap();
    let mut explicit = TransientFilm::new(MATERIAL, 1.0, 17).unwrap();
    explicit.set_wall_slip_length_m(0.0).unwrap();
    for acceleration in [0.4, -0.2, 0.0] {
        let drive = FilmDrive {
            acceleration_m_s2: acceleration,
            wall_velocity_m_s: -0.03,
            surface_shear_pa: 0.7,
        };
        let a = original.step(0.03, drive, LIMITS, || false).unwrap();
        let b = explicit.step(0.03, drive, LIMITS, || false).unwrap();
        assert_eq!(a, b);
        assert_eq!(a.wall_slip_heat_j_m2.to_bits(), 0.0_f64.to_bits());
        for (&u, &v) in original
            .velocities_m_s()
            .iter()
            .zip(explicit.velocities_m_s())
        {
            assert_eq!(u.to_bits(), v.to_bits());
        }
    }
    explicit.set_wall_slip_length_m(0.2).unwrap();
    let state = explicit.velocities_m_s().to_vec();
    let time = explicit.accepted_time_s();
    for length in [-0.1, f64::NAN, f64::INFINITY] {
        assert!(explicit.set_wall_slip_length_m(length).is_err());
        assert_eq!(explicit.wall_slip_length_m(), 0.2);
        assert_eq!(explicit.velocities_m_s(), state);
        assert_eq!(explicit.accepted_time_s(), time);
    }
}

#[test]
fn signed_couette_has_exact_bulk_and_interface_heat_in_every_frame() {
    let (n, b, dt) = (16, 0.3, 0.02);
    let mu = MATERIAL.density_kg_m3 * MATERIAL.kinematic_viscosity_m2_s;
    for rate in [-0.4, 0.4] {
        let mut prior = None;
        for boost in [0.0, 0.7] {
            let wall = -0.2 + boost;
            let initial: Vec<_> = (0..n)
                .map(|i| wall + rate * (b + (i as f64 + 0.5) / n as f64))
                .collect();
            let mut film =
                TransientFilm::from_profile(MATERIAL, 1.0, initial.clone(), 0.0).unwrap();
            film.set_wall_slip_length_m(b).unwrap();
            let drive = FilmDrive {
                wall_velocity_m_s: wall,
                surface_shear_pa: mu * rate,
                ..FilmDrive::default()
            };
            let mut heat = vec![f64::NAN; n];
            let report = film
                .step_with_heat(dt, drive, LIMITS, &mut heat, || false)
                .unwrap();
            for (&u, &v) in film.velocities_m_s().iter().zip(&initial) {
                near(u, v, 2e-15);
            }
            near(
                report.midpoint_wall_fluid_velocity_m_s,
                wall + b * rate,
                2e-15,
            );
            near(report.midpoint_wall_shear_pa, mu * rate, 5e-13);
            near(
                report.final_flow_per_width_m2_s,
                wall + rate * (0.5 + b),
                2e-15,
            );
            near(report.viscous_heat_j_m2, dt * mu * rate * rate, 2e-14);
            near(report.wall_slip_heat_j_m2, dt * mu * rate * rate * b, 2e-14);
            for &q in &heat {
                near(q, dt * mu * rate * rate / n as f64, 2e-14);
            }
            near(heat.iter().sum(), report.viscous_heat_j_m2, 2e-14);
            near(report.wall_work_j_m2, -dt * mu * rate * wall, 2e-14);
            near(report.work_residual_j_m2, 0.0, 2e-12);
            if let Some(first) = prior {
                let first: FilmStep = first;
                near(report.viscous_heat_j_m2, first.viscous_heat_j_m2, 2e-14);
                near(report.wall_slip_heat_j_m2, first.wall_slip_heat_j_m2, 2e-14);
                near(
                    report.wall_work_j_m2 - first.wall_work_j_m2,
                    -dt * mu * rate * boost,
                    2e-14,
                );
                near(
                    report.surface_work_j_m2 - first.surface_work_j_m2,
                    dt * mu * rate * boost,
                    2e-14,
                );
            }
            prior = Some(report);
        }
    }
}

#[test]
fn body_and_opposed_surface_stress_converge_to_independent_slip_flux() {
    let (b, wall, acceleration, traction, dt) = (0.2, -0.1, 0.02, -3.0, 0.03);
    let mu = MATERIAL.density_kg_m3 * MATERIAL.kinematic_viscosity_m2_s;
    let force = MATERIAL.density_kg_m3 * acceleration;
    let tau_wall = force + traction;
    let reference = wall + force / (3.0 * mu) + traction / (2.0 * mu) + b * tau_wall / mu;
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let dy = 1.0 / n as f64;
        // Independent conservation tau(y)=tau_surface+F*(h-y); integrate the
        // admitted wall jump and each native shear segment, not a solver output.
        let mut u = wall + (b + 0.5 * dy) * tau_wall / mu;
        let mut initial = vec![u];
        for i in 1..n {
            u += dy * (traction + force * (1.0 - i as f64 * dy)) / mu;
            initial.push(u);
        }
        let mut film = TransientFilm::from_profile(MATERIAL, 1.0, initial.clone(), 0.0).unwrap();
        film.set_wall_slip_length_m(b).unwrap();
        let report = film
            .step(
                dt,
                FilmDrive {
                    acceleration_m_s2: acceleration,
                    wall_velocity_m_s: wall,
                    surface_shear_pa: traction,
                },
                LIMITS,
                || false,
            )
            .unwrap();
        for (&u, &v) in film.velocities_m_s().iter().zip(&initial) {
            near(u, v, 2e-14);
        }
        near(
            report.wall_slip_heat_j_m2,
            dt * b * tau_wall * tau_wall / mu,
            2e-13,
        );
        near(report.work_residual_j_m2, 0.0, 2e-12);
        errors.push((report.final_flow_per_width_m2_s - reference).abs());
    }
    assert!(errors[0] / errors[1] > 3.99, "{errors:?}");
    assert!(errors[1] / errors[2] > 3.99, "{errors:?}");
    assert!(errors[2] < 1.1e-4, "{errors:?}");
}

#[test]
fn two_cell_midpoint_matches_closed_form_and_late_refusal_keeps_the_boundary() {
    let material = FilmMaterial {
        density_kg_m3: 2.0,
        kinematic_viscosity_m2_s: 0.3,
        surface_tension_n_m: 0.0,
    };
    let (b, dt, dy) = (0.2, 0.2, 0.5);
    let drive = FilmDrive {
        acceleration_m_s2: -0.4,
        wall_velocity_m_s: 0.3,
        surface_shear_pa: 0.1,
    };
    let initial = [0.2, -0.1];
    let make = || {
        let mut film = TransientFilm::from_profile(material, 1.0, initial.to_vec(), 0.0).unwrap();
        film.set_wall_slip_length_m(b).unwrap();
        film
    };
    // Direct 2x2 midpoint equations, independently solved by their determinant.
    let (mass, mu) = (
        material.density_kg_m3 * dy,
        material.density_kg_m3 * material.kinematic_viscosity_m2_s,
    );
    let (k, kw) = (mu / dy, mu / (0.5 * dy + b));
    let (aa, ab, bb) = (mass / dt + 0.5 * (k + kw), -0.5 * k, mass / dt + 0.5 * k);
    let r0 = mass * drive.acceleration_m_s2 + k * (initial[1] - initial[0])
        - kw * (initial[0] - drive.wall_velocity_m_s);
    let r1 =
        mass * drive.acceleration_m_s2 + drive.surface_shear_pa - k * (initial[1] - initial[0]);
    let determinant = aa * bb - ab * ab;
    let exact = [
        initial[0] + (bb * r0 - ab * r1) / determinant,
        initial[1] + (aa * r1 - ab * r0) / determinant,
    ];
    let mut reference = make();
    let mut heat = [0.0; 2];
    let mut polls = 0;
    let expected = reference
        .step_with_heat(dt, drive, LIMITS, &mut heat, || {
            polls += 1;
            false
        })
        .unwrap();
    for (&u, &v) in reference.velocities_m_s().iter().zip(&exact) {
        near(u, v, 2e-15);
    }
    near(expected.work_residual_j_m2, 0.0, 2e-15);
    let mut film = make();
    let mut output = [-7.0; 2];
    let mut seen = 0;
    assert_eq!(
        film.step_with_heat(dt, drive, LIMITS, &mut output, || {
            seen += 1;
            seen == polls
        }),
        Err(TransientFilmError::Cancelled)
    );
    assert_eq!(film.velocities_m_s(), initial);
    assert_eq!(film.accepted_time_s(), 0.0);
    assert_eq!(film.wall_slip_length_m(), b);
    assert_eq!(output, [-7.0; 2]);
    assert_eq!(
        film.step_with_heat(dt, drive, LIMITS, &mut output, || false)
            .unwrap(),
        expected
    );
    assert_eq!(output, heat);
}

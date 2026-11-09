//! Original local material, independent steady stress, conservative heat and
//! transactional publication. Synthetic laws do not calibrate a real fluid.

use super::*;

fn temperatures(cells: usize) -> Vec<f64> {
    (0..=cells)
        .map(|i| 290.0 + 50.0 * i as f64 / cells as f64)
        .collect()
}

// No PreparedCarreau inverse, slit profile or numerical film owner enters this
// reference. Monotone bisection evaluates the original complete thermal law.
fn independent_rate(source: &ThermalCarreauYasuda, temperature: f64, stress: f64) -> f64 {
    if stress == 0.0 {
        return 0.0;
    }
    let target = stress.abs();
    let evaluate = |rate| {
        source
            .evaluate(temperature, rate)
            .unwrap()
            .dynamic_viscosity_pa_s()
            * rate
    };
    let (mut low, mut high) = (0.0, SPEC.maximum_laminar_shear_rate_per_s);
    assert!(evaluate(high) >= target);
    for _ in 0..80 {
        let midpoint = 0.5 * (low + high);
        if midpoint == low || midpoint == high {
            break;
        }
        if evaluate(midpoint) < target {
            low = midpoint;
        } else {
            high = midpoint;
        }
    }
    (0.5 * (low + high)).copysign(stress)
}

// For h=rho=1, steady momentum is tau(y)=tau_surface+g*(1-y).
// Simpson integrates gamma, (1-y)*gamma and tau*gamma independently of the
// production finite volumes, giving surface speed, volume flux and dissipation.
fn continuous_reference(source: &ThermalCarreauYasuda, drive: FilmDrive, n: usize) -> [f64; 3] {
    let mut integrals = [0.0; 3];
    for j in 0..=n {
        let y = j as f64 / n as f64;
        let stress = drive.surface_shear_pa + drive.acceleration_m_s2 * (1.0 - y);
        let rate = independent_rate(source, 290.0 + 50.0 * y, stress);
        let weight = if j == 0 || j == n {
            1.0
        } else if j % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for (sum, value) in integrals
            .iter_mut()
            .zip([rate, (1.0 - y) * rate, stress * rate])
        {
            *sum += weight * value / (3.0 * n as f64);
        }
    }
    [
        drive.wall_velocity_m_s + integrals[0],
        drive.wall_velocity_m_s + integrals[1],
        integrals[2],
    ]
}

#[test]
fn linear_temperature_steady_stress_converges_to_independent_continuous_flow_and_heat() {
    let source = law(0.45);
    let drive = FilmDrive {
        acceleration_m_s2: 3.6,
        wall_velocity_m_s: 0.03,
        surface_shear_pa: -1.2,
    };
    let coarse_reference = continuous_reference(&source, drive, 1024);
    let reference = continuous_reference(&source, drive, 2048);
    for (a, b) in coarse_reference.into_iter().zip(reference) {
        near(a, b, 2e-11);
    }
    near(
        drive.acceleration_m_s2 * reference[1] + drive.surface_shear_pa * reference[0]
            - (drive.surface_shear_pa + drive.acceleration_m_s2) * drive.wall_velocity_m_s,
        reference[2],
        2e-13,
    );
    let mut previous_error = [f64::INFINITY; 3];
    for cells in [16, 32, 64] {
        let temperature = temperatures(cells);
        let spacing = 1.0 / cells as f64;
        let stress: Vec<_> = (0..=cells)
            .map(|i| drive.surface_shear_pa + drive.acceleration_m_s2 * (1.0 - i as f64 * spacing))
            .collect();
        let rates: Vec<_> = temperature
            .iter()
            .zip(&stress)
            .map(|(&t, &tau)| independent_rate(&source, t, tau))
            .collect();
        assert!(rates[0] > 0.0 && rates[cells] < 0.0);
        let mut velocity = drive.wall_velocity_m_s + 0.5 * spacing * rates[0];
        let mut initial = vec![velocity];
        for &rate in &rates[1..cells] {
            velocity += spacing * rate;
            initial.push(velocity);
        }
        let mut film =
            TransientCarreauFilm::from_profile(source, 300.0, SPEC, initial.clone(), 0.0).unwrap();
        film.set_face_temperatures(&temperature, || false).unwrap();
        let mut heat = vec![f64::NAN; cells];
        let dt = 0.05;
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
        for (&actual, &expected) in film.velocities_m_s().iter().zip(&initial) {
            near(actual, expected, 2e-12);
        }
        for i in 0..cells {
            let expected =
                0.5 * dt * spacing * (stress[i] * rates[i] + stress[i + 1] * rates[i + 1]);
            near(heat[i], expected, 2e-13);
            assert!(heat[i] >= 0.0);
        }
        near(heat.iter().sum(), report.flow.viscous_heat_j_m2, 3e-13);
        near(report.flow.work_residual_j_m2, 0.0, 3e-13);
        let actual = [
            report.flow.midpoint_surface_velocity_m_s,
            report.flow.final_flow_per_width_m2_s,
            report.flow.viscous_heat_j_m2 / dt,
        ];
        let error = std::array::from_fn(|i| (actual[i] - reference[i]).abs());
        for i in 0..3 {
            assert!(
                error[i] < 0.27 * previous_error[i],
                "cells={cells}, errors={error:?}, previous={previous_error:?}"
            );
        }
        previous_error = error;
    }
    assert!(previous_error[0] < 9e-5);
    assert!(previous_error[1] < 1.5e-4);
    assert!(previous_error[2] < 7e-4);
}

#[test]
fn uniform_face_temperatures_preserve_original_step_bits_and_restore_the_compact_source() {
    let cells = 16;
    let mut original = TransientCarreauFilm::new(law(0.45), 330.0, SPEC, cells).unwrap();
    let mut local = TransientCarreauFilm::new(law(0.45), 300.0, SPEC, cells).unwrap();
    local
        .set_face_temperatures(&vec![330.0; cells + 1], || false)
        .unwrap();
    assert!(local.is_isothermal());
    assert!(local.material_at_face(cells + 1).is_none());
    let drive = FilmDrive {
        acceleration_m_s2: 0.7,
        surface_shear_pa: -0.1,
        ..FilmDrive::default()
    };
    let mut heat = vec![0.0; cells];
    for _ in 0..3 {
        let a = original
            .step(
                0.03,
                drive,
                LIMITS,
                TransientCarreauControls::default(),
                || false,
            )
            .unwrap();
        let b = local
            .step_with_heat(
                0.03,
                drive,
                LIMITS,
                TransientCarreauControls::default(),
                &mut heat,
                || false,
            )
            .unwrap();
        assert_eq!(a, b);
        assert_eq!(original.velocities_m_s(), local.velocities_m_s());
        near(heat.iter().sum(), b.flow.viscous_heat_j_m2, 1e-13);
    }
    let before = local.velocities_m_s().to_vec();
    let time = local.accepted_time_s();
    local
        .set_face_temperatures(&temperatures(cells), || false)
        .unwrap();
    assert!(!local.is_isothermal());
    local.set_temperature(330.0, || false).unwrap();
    assert!(local.is_isothermal());
    for i in 0..=cells {
        assert_eq!(local.material_at_face(i).unwrap(), original.material());
    }
    assert_eq!(local.velocities_m_s(), before);
    assert_eq!(local.accepted_time_s(), time);
}

#[test]
fn spatial_transient_heat_uses_each_original_face_law_and_is_frame_invariant() {
    let (cells, boost, dt) = (24, 1.25, 0.03);
    let source = law(0.45);
    let temperature = temperatures(cells);
    let mut film = TransientCarreauFilm::new(source, 300.0, SPEC, cells).unwrap();
    let mut moving =
        TransientCarreauFilm::from_profile(source, 300.0, SPEC, vec![boost; cells], 0.0).unwrap();
    for owner in [&mut film, &mut moving] {
        owner.set_face_temperatures(&temperature, || false).unwrap();
    }
    let mut heat = vec![0.0; cells];
    let mut boosted_heat = vec![0.0; cells];
    for (acceleration, traction) in [(0.8, 0.1), (-0.5, -0.2), (0.0, 0.3)] {
        let initial = film.velocities_m_s().to_vec();
        let drive = FilmDrive {
            acceleration_m_s2: acceleration,
            wall_velocity_m_s: 0.03,
            surface_shear_pa: traction,
        };
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
        let shifted = moving
            .step_with_heat(
                dt,
                FilmDrive {
                    wall_velocity_m_s: drive.wall_velocity_m_s + boost,
                    ..drive
                },
                LIMITS,
                TransientCarreauControls::default(),
                &mut boosted_heat,
                || false,
            )
            .unwrap();
        let midpoint: Vec<_> = initial
            .iter()
            .zip(film.velocities_m_s())
            .map(|(&a, &b)| 0.5 * (a + b))
            .collect();
        let spacing = film.spacing_m();
        let mut segment = Vec::with_capacity(cells + 1);
        for face in 0..=cells {
            let (difference, distance) = if face == 0 {
                (midpoint[0] - drive.wall_velocity_m_s, 0.5 * spacing)
            } else if face == cells {
                (
                    report.flow.midpoint_surface_velocity_m_s - midpoint[cells - 1],
                    0.5 * spacing,
                )
            } else {
                (midpoint[face] - midpoint[face - 1], spacing)
            };
            let rate = difference / distance;
            let mu = source
                .evaluate(temperature[face], rate.abs())
                .unwrap()
                .dynamic_viscosity_pa_s();
            segment.push(dt * mu * rate * difference);
        }
        for i in 0..cells {
            let expected = segment[i] * if i == 0 { 1.0 } else { 0.5 }
                + segment[i + 1] * if i + 1 == cells { 1.0 } else { 0.5 };
            near(heat[i], expected, 3e-13);
            near(heat[i], boosted_heat[i], 3e-13);
            assert!(heat[i] >= 0.0);
            near(
                moving.velocities_m_s()[i] - film.velocities_m_s()[i],
                boost,
                3e-12,
            );
        }
        near(heat.iter().sum(), report.flow.viscous_heat_j_m2, 3e-13);
        near(
            report.flow.viscous_heat_j_m2,
            shifted.flow.viscous_heat_j_m2,
            3e-13,
        );
        near(report.flow.work_residual_j_m2, 0.0, 2e-12);
    }
}

#[test]
fn local_temperature_refusal_and_late_cancellation_preserve_every_accepted_face() {
    let cells = 16;
    let make = || TransientCarreauFilm::new(law(0.45), 300.0, SPEC, cells).unwrap();
    let mut film = make();
    let initial = film.velocities_m_s().to_vec();
    let initial_material = *film.material();
    let initial_scale = film.zero_shear_diffusion_time_s();
    let field = temperatures(cells);
    let mut invalid = field.clone();
    invalid[cells] = 361.0;
    for rejected in [&invalid[..], &field[..cells]] {
        assert!(film.set_face_temperatures(rejected, || false).is_err());
        assert!(film.is_isothermal());
        for face in 0..=cells {
            assert_eq!(film.material_at_face(face).unwrap(), &initial_material);
        }
    }
    let mut reference = make();
    let mut polls = 0;
    reference
        .set_face_temperatures(&field, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        film.set_face_temperatures(&field, || {
            calls += 1;
            calls == polls
        }),
        Err(TransientCarreauError::Flow(TransientFilmError::Cancelled))
    ));
    assert!(film.is_isothermal());
    assert_eq!(*film.material(), initial_material);
    assert_eq!(film.zero_shear_diffusion_time_s(), initial_scale);
    assert_eq!(film.velocities_m_s(), initial);
    assert_eq!(film.accepted_time_s(), 0.0);
    film.set_face_temperatures(&field, || false).unwrap();
    for (face, &t) in field.iter().enumerate() {
        let actual = film.material_at_face(face).unwrap();
        assert_eq!(actual, reference.material_at_face(face).unwrap());
        assert_eq!(actual.temperature_k(), t);
        assert_eq!(actual.law(), &law(0.45));
    }
}

#[test]
fn heat_output_and_history_publish_together_after_all_refusal_and_cancellation_checks() {
    let cells = 16;
    let make = || {
        let mut film = TransientCarreauFilm::new(law(0.45), 300.0, SPEC, cells).unwrap();
        film.set_face_temperatures(&temperatures(cells), || false)
            .unwrap();
        film
    };
    let drive = FilmDrive {
        acceleration_m_s2: 0.8,
        surface_shear_pa: -0.1,
        ..FilmDrive::default()
    };
    let controls = TransientCarreauControls::default();
    let mut reference = make();
    let mut expected_heat = vec![0.0; cells];
    let mut polls = 0;
    let expected = reference
        .step_with_heat(0.05, drive, LIMITS, controls, &mut expected_heat, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut victim = make();
    let initial = victim.velocities_m_s().to_vec();
    let mut short = vec![7.0; cells - 1];
    assert!(
        victim
            .step_with_heat(0.05, drive, LIMITS, controls, &mut short, || false)
            .is_err()
    );
    assert!(short.iter().all(|&q| q == 7.0));
    let mut heat = vec![7.0; cells];
    assert!(
        victim
            .step_with_heat(
                0.05,
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
    let mut calls = 0;
    assert!(matches!(
        victim.step_with_heat(0.05, drive, LIMITS, controls, &mut heat, || {
            calls += 1;
            calls == polls
        }),
        Err(TransientCarreauError::Flow(TransientFilmError::Cancelled))
    ));
    assert!(heat.iter().all(|&q| q == 7.0));
    assert_eq!(victim.velocities_m_s(), initial);
    assert_eq!(victim.accepted_time_s(), 0.0);
    let accepted = victim
        .step_with_heat(0.05, drive, LIMITS, controls, &mut heat, || false)
        .unwrap();
    assert_eq!(accepted, expected);
    assert_eq!(heat, expected_heat);
    assert_eq!(victim.velocities_m_s(), reference.velocities_m_s());
}

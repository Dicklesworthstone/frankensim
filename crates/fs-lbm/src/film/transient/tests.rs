//! G1 analytic mixed-boundary diffusion; G0/G3 work, forcing and frame laws;
//! G4 transaction refusal. These are numerical model tests, not a measured
//! fluid calibration, turbulent-transition criterion or free-surface-wave test.

use super::*;

#[path = "tests/slip.rs"]
mod slip;

const MATERIAL: FilmMaterial = FilmMaterial {
    density_kg_m3: 1000.0,
    kinematic_viscosity_m2_s: 0.03,
    surface_tension_n_m: 0.072,
};
const LIMITS: FilmStepLimits = FilmStepLimits {
    maximum_speed_m_s: 10.0,
    maximum_velocity_residual_m_s: 1e-10,
};

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}

fn startup_velocity(y: f64, time: f64, acceleration: f64) -> f64 {
    // Independent continuous separation of variables: u(0)=0, u_y(1)=0,
    // initially at rest; eigenvalues k_j=(j+1/2)*pi.
    let nu = MATERIAL.kinematic_viscosity_m2_s;
    let mut value = acceleration * (y - 0.5 * y * y) / nu;
    for j in 0..100 {
        let k = (j as f64 + 0.5) * core::f64::consts::PI;
        value -= 2.0 * acceleration / (nu * k.powi(3)) * (k * y).sin() * (-nu * k * k * time).exp();
    }
    value
}

#[test]
fn startup_resolves_viscous_memory_and_converges_to_the_continuous_solution() {
    let mut errors = Vec::new();
    for n in [16, 32, 64] {
        let mut film = TransientFilm::new(MATERIAL, 1.0, n).unwrap();
        let dt = 0.01;
        for _ in 0..100 {
            let step = film
                .step(
                    dt,
                    FilmDrive {
                        acceleration_m_s2: 0.1,
                        ..FilmDrive::default()
                    },
                    LIMITS,
                    || false,
                )
                .unwrap();
            near(step.work_residual_j_m2, 0.0, 1e-11);
            near(step.momentum_residual_kg_m_s, 0.0, 1e-11);
            assert!(step.viscous_heat_j_m2 >= 0.0);
        }
        let error = film
            .velocities_m_s()
            .iter()
            .enumerate()
            .map(|(i, &u)| {
                let y = (i as f64 + 0.5) / n as f64;
                (u - startup_velocity(y, 1.0, 0.1)).powi(2)
            })
            .sum::<f64>()
            .sqrt()
            / (n as f64).sqrt();
        // Startup remains far from the instantaneous steady parabolic value.
        assert!(*film.velocities_m_s().last().unwrap() < 0.11);
        errors.push(error);
    }
    assert!(errors[0] / errors[1] > 3.7, "{errors:?}");
    assert!(errors[1] / errors[2] > 3.7, "{errors:?}");
    assert!(errors[2] < 3e-5, "{errors:?}");
}

#[test]
fn midpoint_has_second_order_time_convergence_for_the_shear_free_eigenmode() {
    let n = 256;
    let k = core::f64::consts::FRAC_PI_2;
    let profile: Vec<_> = (0..n)
        .map(|i| (k * (i as f64 + 0.5) / n as f64).sin())
        .collect();
    let expected = (-MATERIAL.kinematic_viscosity_m2_s * k * k * 4.0).exp();
    let mut errors = Vec::new();
    for steps in [2, 4, 8] {
        let mut film = TransientFilm::from_profile(MATERIAL, 1.0, profile.clone(), 0.0).unwrap();
        for _ in 0..steps {
            let receipt = film
                .step(4.0 / steps as f64, FilmDrive::default(), LIMITS, || false)
                .unwrap();
            near(
                receipt.kinetic_energy_change_j_m2 + receipt.viscous_heat_j_m2,
                0.0,
                2e-9,
            );
            assert!(receipt.kinetic_energy_change_j_m2 < 0.0);
        }
        errors.push((film.velocities_m_s()[n - 1] / profile[n - 1] - expected).abs());
    }
    assert!(errors[0] / errors[1] > 3.9, "{errors:?}");
    assert!(errors[1] / errors[2] > 3.9, "{errors:?}");
}

#[test]
fn signed_wall_and_surface_work_close_through_reversal_and_common_frame_boost() {
    let n = 40;
    let boost = 0.7;
    let initial: Vec<_> = (0..n).map(|i| 0.04 * (i as f64 + 0.5) / n as f64).collect();
    let boosted = initial.iter().map(|u| u + boost).collect();
    let mut base = TransientFilm::from_profile(MATERIAL, 1.0, initial, 0.0).unwrap();
    let mut moving = TransientFilm::from_profile(MATERIAL, 1.0, boosted, 0.0).unwrap();
    for acceleration in [0.2, -0.3, 0.0] {
        let drive = FilmDrive {
            acceleration_m_s2: acceleration,
            wall_velocity_m_s: 0.01,
            surface_shear_pa: -0.2,
        };
        let a = base.step(0.05, drive, LIMITS, || false).unwrap();
        let b = moving
            .step(
                0.05,
                FilmDrive {
                    wall_velocity_m_s: drive.wall_velocity_m_s + boost,
                    ..drive
                },
                LIMITS,
                || false,
            )
            .unwrap();
        for (&u, &v) in base.velocities_m_s().iter().zip(moving.velocities_m_s()) {
            near(v - u, boost, 2e-15);
        }
        near(a.viscous_heat_j_m2, b.viscous_heat_j_m2, 2e-13);
        near(a.midpoint_wall_shear_pa, b.midpoint_wall_shear_pa, 2e-12);
        near(a.work_residual_j_m2, 0.0, 2e-12);
        near(b.work_residual_j_m2, 0.0, 2e-11);
        near(
            b.kinetic_energy_change_j_m2 - a.kinetic_energy_change_j_m2,
            boost * a.momentum_change_kg_m_s,
            2e-11,
        );
        assert!(a.surface_work_j_m2 < 0.0);
        assert!(a.viscous_heat_j_m2 > 0.0);
    }
}

#[test]
fn steady_nusselt_and_surface_traction_limit_converges_with_actual_boundary_half_cells() {
    let drive = FilmDrive {
        acceleration_m_s2: -0.02,
        wall_velocity_m_s: 0.05,
        surface_shear_pa: 3.0,
    };
    let steady = ShearedFilm::new(
        MATERIAL,
        1.0,
        drive.acceleration_m_s2,
        drive.surface_shear_pa,
    )
    .unwrap();
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let profile = (0..n)
            .map(|i| {
                drive.wall_velocity_m_s + steady.velocity_m_s((i as f64 + 0.5) / n as f64).unwrap()
            })
            .collect();
        let mut film = TransientFilm::from_profile(MATERIAL, 1.0, profile, 0.0).unwrap();
        let mut last = None;
        for _ in 0..200 {
            last = Some(film.step(1.0, drive, LIMITS, || false).unwrap());
        }
        let step = last.unwrap();
        near(step.work_residual_j_m2, 0.0, 2e-11);
        near(step.midpoint_wall_shear_pa, steady.wall_shear_pa(), 2e-6);
        let analytic_flux = steady.flow_per_width_m2_s() + drive.wall_velocity_m_s;
        errors.push((step.final_flow_per_width_m2_s - analytic_flux).abs());
    }
    assert!(errors[0] / errors[1] > 3.9, "{errors:?}");
    assert!(errors[1] / errors[2] > 3.9, "{errors:?}");
}

#[test]
fn late_cancellation_and_speed_refusal_preserve_history_and_exact_retry() {
    let mut film = TransientFilm::new(MATERIAL, 1.0, 257).unwrap();
    let drive = FilmDrive {
        acceleration_m_s2: 0.3,
        ..FilmDrive::default()
    };
    film.step(0.01, drive, LIMITS, || false).unwrap();
    let saved = film.velocities_m_s().to_vec();
    let time = film.accepted_time_s();
    let mut reference = TransientFilm::from_profile(MATERIAL, 1.0, saved.clone(), time).unwrap();
    let expected = reference.step(0.01, drive, LIMITS, || false).unwrap();
    // Measure the real warmed trajectory's final poll rather than guessing a
    // count tied to internal linear-solver loops.
    let mut probe = TransientFilm::from_profile(MATERIAL, 1.0, saved.clone(), time).unwrap();
    probe.step(0.01, drive, LIMITS, || false).unwrap();
    let mut calls = 0;
    probe
        .step(0.01, drive, LIMITS, || {
            calls += 1;
            false
        })
        .unwrap();
    let mut seen = 0;
    assert_eq!(
        film.step(0.01, drive, LIMITS, || {
            seen += 1;
            seen == calls
        }),
        Err(TransientFilmError::Cancelled)
    );
    assert_eq!(film.velocities_m_s(), saved);
    assert_eq!(film.accepted_time_s(), time);
    let too_slow = FilmStepLimits {
        maximum_speed_m_s: 1e-6,
        ..LIMITS
    };
    assert_eq!(
        film.step(0.01, drive, too_slow, || false),
        Err(TransientFilmError::SpeedLimit)
    );
    assert_eq!(film.velocities_m_s(), saved);
    let retry = film.step(0.01, drive, LIMITS, || false).unwrap();
    assert_eq!(retry, expected);
    assert_eq!(film.velocities_m_s(), reference.velocities_m_s());
}

#[test]
fn initial_and_drive_domains_refuse_without_clipping() {
    assert!(TransientFilm::new(MATERIAL, 0.0, 8).is_err());
    assert!(TransientFilm::new(MATERIAL, 1.0, 1).is_err());
    assert!(TransientFilm::new(MATERIAL, 1.0, MAX_FILM_CELLS + 1).is_err());
    assert!(TransientFilm::from_profile(MATERIAL, 1.0, vec![0.0, f64::NAN], 0.0).is_err());
    let mut film = TransientFilm::new(MATERIAL, 1.0, 8).unwrap();
    assert!(
        film.step(0.0, FilmDrive::default(), LIMITS, || false)
            .is_err()
    );
    assert!(
        film.step(
            0.1,
            FilmDrive {
                surface_shear_pa: f64::NAN,
                ..FilmDrive::default()
            },
            LIMITS,
            || false
        )
        .is_err()
    );
    assert_eq!(film.accepted_time_s(), 0.0);
    assert!(film.velocities_m_s().iter().all(|&u| u == 0.0));
    let rest = film
        .step(0.1, FilmDrive::default(), LIMITS, || false)
        .unwrap();
    assert_eq!(rest.viscous_heat_j_m2, 0.0);
    assert_eq!(rest.final_flow_per_width_m2_s, 0.0);
}

#[test]
fn resolved_cell_heat_recovers_exact_signed_couette_heating_and_native_step() {
    let n = 16;
    let dt = 0.02;
    let mu = MATERIAL.density_kg_m3 * MATERIAL.kinematic_viscosity_m2_s;
    for rate in [-0.4, 0.4] {
        let drive = FilmDrive {
            wall_velocity_m_s: -0.2,
            surface_shear_pa: mu * rate,
            ..FilmDrive::default()
        };
        let profile: Vec<_> = (0..n)
            .map(|i| drive.wall_velocity_m_s + rate * (i as f64 + 0.5) / n as f64)
            .collect();
        let mut plain = TransientFilm::from_profile(MATERIAL, 1.0, profile.clone(), 0.0).unwrap();
        let mut thermal = TransientFilm::from_profile(MATERIAL, 1.0, profile, 0.0).unwrap();
        let reference = plain.step(dt, drive, LIMITS, || false).unwrap();
        let mut heat = vec![f64::NAN; n];
        let report = thermal
            .step_with_heat(dt, drive, LIMITS, &mut heat, || false)
            .unwrap();
        assert_eq!(report, reference);
        assert_eq!(thermal.velocities_m_s(), plain.velocities_m_s());
        for value in &heat {
            near(*value, dt * mu * rate * rate / n as f64, 2e-16);
        }
        near(heat.iter().sum(), report.viscous_heat_j_m2, 2e-16);
        near(report.kinetic_energy_change_j_m2, 0.0, 1e-13);
        near(report.work_residual_j_m2, 0.0, 1e-13);
    }
}

#[test]
fn heat_output_and_flow_refuse_and_retry_as_one_transaction() {
    let n = 257;
    let drive = FilmDrive {
        acceleration_m_s2: 0.3,
        wall_velocity_m_s: 0.01,
        surface_shear_pa: -0.02,
    };
    let mut film = TransientFilm::new(MATERIAL, 1.0, n).unwrap();
    let mut reference = TransientFilm::new(MATERIAL, 1.0, n).unwrap();
    film.step(0.01, drive, LIMITS, || false).unwrap();
    reference.step(0.01, drive, LIMITS, || false).unwrap();
    let original = film.velocities_m_s().to_vec();
    let original_time = film.accepted_time_s();
    let mut expected_heat = vec![0.0; n];
    let mut polls = 0;
    let expected = reference
        .step_with_heat(0.01, drive, LIMITS, &mut expected_heat, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut heat = vec![-17.0; n];
    assert!(
        film.step_with_heat(0.01, drive, LIMITS, &mut heat[..n - 1], || false)
            .is_err()
    );
    let mut seen = 0;
    assert_eq!(
        film.step_with_heat(0.01, drive, LIMITS, &mut heat, || {
            seen += 1;
            seen == polls
        }),
        Err(TransientFilmError::Cancelled)
    );
    assert_eq!(film.velocities_m_s(), original);
    assert_eq!(film.accepted_time_s(), original_time);
    assert!(heat.iter().all(|&q| q == -17.0));
    let accepted = film
        .step_with_heat(0.01, drive, LIMITS, &mut heat, || false)
        .unwrap();
    assert_eq!(accepted, expected);
    assert_eq!(heat, expected_heat);
    assert_eq!(film.velocities_m_s(), reference.velocities_m_s());
}

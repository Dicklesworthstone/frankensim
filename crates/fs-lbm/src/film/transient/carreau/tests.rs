//! G1 constitutive/steady/temporal limits; G3 signed work and frame covariance;
//! G4 source/work/cancellation refusal. No measured-fluid or transition claim.

use super::*;
use fs_material::fluid_transport::rheology::{CarreauYasuda, ViscosityShiftMode};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

#[path = "tests/thermal.rs"]
mod thermal;

#[path = "tests/slip.rs"]
mod slip;

const SPEC: TransientCarreauFilmSpec = TransientCarreauFilmSpec {
    density_kg_m3: 1.0,
    thickness_m: 1.0,
    maximum_laminar_shear_rate_per_s: 20.0,
};
const LIMITS: FilmStepLimits = FilmStepLimits {
    maximum_speed_m_s: 20.0,
    maximum_velocity_residual_m_s: 1e-12,
};

fn law(exponent: f64) -> ThermalCarreauYasuda {
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(2.0, 0.25, 2.0, exponent, 2.0, [0.0, 1e5]).unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 12_000.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        ViscosityShiftMode::VerticalAndRate,
    )
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn prepared_material_preserves_original_temperature_rate_and_positive_tangent() {
    let source = law(0.45);
    for temperature in [290.0, 300.0, 330.0] {
        let prepared = PreparedCarreau::new(source, temperature, 20.0).unwrap();
        assert_eq!(prepared.law(), &source);
        assert_eq!(prepared.temperature_k(), temperature);
        for rate in [0.0, 0.1, 3.0, 12.0, 20.0] {
            let state = source.evaluate(temperature, rate).unwrap();
            let [mu, tangent] = prepared.viscosity_and_tangent(rate).unwrap();
            near(mu, state.dynamic_viscosity_pa_s(), 2e-14 * mu);
            near(tangent, state.stress_tangent_pa_s(), 2e-14 * tangent);
            for sign in [-1.0, 1.0] {
                let (inverse, _) = prepared
                    .signed_shear_rate_at_stress(sign * mu * rate, CarreauPipeControls::default())
                    .unwrap();
                near(inverse, sign * rate, 2e-10 * rate.max(1e-10));
            }
        }
        assert!(prepared.viscosity_and_tangent(20.1).is_err());
        assert!(
            prepared
                .signed_shear_rate_at_stress(f64::NAN, CarreauPipeControls::default())
                .is_err()
        );
    }
}

#[test]
fn constant_source_recovers_the_original_transient_newtonian_solver() {
    for temperature in [300.0, 330.0] {
        let source = law(1.0);
        let mu = source
            .evaluate(temperature, 0.0)
            .unwrap()
            .dynamic_viscosity_pa_s();
        let mut nonlinear = TransientCarreauFilm::new(source, temperature, SPEC, 24).unwrap();
        let mut linear = TransientFilm::new(
            FilmMaterial {
                density_kg_m3: SPEC.density_kg_m3,
                kinematic_viscosity_m2_s: mu / SPEC.density_kg_m3,
                surface_tension_n_m: 0.0,
            },
            SPEC.thickness_m,
            24,
        )
        .unwrap();
        for acceleration in [0.4, -0.7, 0.0] {
            for _ in 0..20 {
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
                let b = linear.step(0.02, drive, LIMITS, || false).unwrap();
                for (&u, &v) in nonlinear
                    .velocities_m_s()
                    .iter()
                    .zip(linear.velocities_m_s())
                {
                    near(u, v, 2e-12);
                }
                near(a.flow.viscous_heat_j_m2, b.viscous_heat_j_m2, 1e-12);
                near(a.flow.work_residual_j_m2, 0.0, 2e-12);
                assert!(a.newton_iterations <= 1);
            }
        }
    }
}

// Independent reference formula at T_ref, with bisection rather than the
// production constitutive inverse. Integrate its stress-selected rate in y.
fn reference_rate(stress: f64) -> f64 {
    let mut low = 0.0;
    let mut high = 20.0;
    for _ in 0..80 {
        let rate = 0.5 * (low + high);
        let mu = 0.25 + 1.75 * (1.0 + (2.0_f64 * rate).powi(2)).powf(-0.25);
        if mu * rate < stress.abs() {
            low = rate;
        } else {
            high = rate;
        }
    }
    (0.5 * (low + high)).copysign(stress)
}

fn reference_velocity(y: f64) -> f64 {
    let panels = 1024;
    let dy = y / panels as f64;
    let mut integral = 0.0;
    for i in 0..=panels {
        let stress = 0.15 + 0.8 * (1.0 - i as f64 * dy);
        let weight = if i == 0 || i == panels {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        integral += weight * reference_rate(stress);
    }
    dy * integral / 3.0
}

#[test]
fn nonlinear_startup_converges_to_independently_integrated_steady_flow() {
    let mut errors = Vec::new();
    for cells in [8, 16, 32] {
        let mut film = TransientCarreauFilm::new(law(0.5), 300.0, SPEC, cells).unwrap();
        let drive = FilmDrive {
            acceleration_m_s2: 0.8,
            surface_shear_pa: 0.15,
            ..FilmDrive::default()
        };
        let first = film
            .step(
                0.025,
                drive,
                LIMITS,
                TransientCarreauControls::default(),
                || false,
            )
            .unwrap();
        assert!(first.flow.final_flow_per_width_m2_s < 0.03);
        for _ in 1..1200 {
            film.step(
                0.025,
                drive,
                LIMITS,
                TransientCarreauControls::default(),
                || false,
            )
            .unwrap();
        }
        let error = film
            .velocities_m_s()
            .iter()
            .enumerate()
            .map(|(i, &u)| {
                let y = (i as f64 + 0.5) / cells as f64;
                (u - reference_velocity(y)).powi(2)
            })
            .sum::<f64>()
            .sqrt()
            / (cells as f64).sqrt();
        errors.push(error);
    }
    assert!(errors[0] / errors[1] > 3.7, "{errors:?}");
    assert!(errors[1] / errors[2] > 3.7, "{errors:?}");
    assert!(errors[2] < 2e-4, "{errors:?}");
}

#[test]
fn nonlinear_reversal_and_moving_wall_preserve_work_and_frame_covariance() {
    let cells = 20;
    let boost = 0.6;
    let initial: Vec<_> = (0..cells)
        .map(|i| 0.02 * (i as f64 + 0.5) / cells as f64)
        .collect();
    let mut base =
        TransientCarreauFilm::from_profile(law(0.4), 300.0, SPEC, initial.clone(), 0.0).unwrap();
    let mut moving = TransientCarreauFilm::from_profile(
        law(0.4),
        300.0,
        SPEC,
        initial.iter().map(|u| u + boost).collect(),
        0.0,
    )
    .unwrap();
    let mut final_q = 0.0;
    for acceleration in [0.7, -0.9] {
        for _ in 0..60 {
            let drive = FilmDrive {
                acceleration_m_s2: acceleration,
                wall_velocity_m_s: 0.01,
                surface_shear_pa: -0.03,
            };
            let a = base
                .step(
                    0.05,
                    drive,
                    LIMITS,
                    TransientCarreauControls::default(),
                    || false,
                )
                .unwrap();
            let b = moving
                .step(
                    0.05,
                    FilmDrive {
                        wall_velocity_m_s: drive.wall_velocity_m_s + boost,
                        ..drive
                    },
                    LIMITS,
                    TransientCarreauControls::default(),
                    || false,
                )
                .unwrap();
            for (&u, &v) in base.velocities_m_s().iter().zip(moving.velocities_m_s()) {
                near(v - u, boost, 2e-11);
            }
            near(a.flow.viscous_heat_j_m2, b.flow.viscous_heat_j_m2, 2e-11);
            near(a.flow.work_residual_j_m2, 0.0, 2e-11);
            near(b.flow.work_residual_j_m2, 0.0, 2e-11);
            near(a.flow.momentum_residual_kg_m_s, 0.0, 2e-11);
            assert!(a.flow.viscous_heat_j_m2 >= 0.0);
            assert!(a.surface_stress_residual_pa.abs() < 1e-12);
            final_q = a.flow.final_flow_per_width_m2_s;
        }
    }
    assert!(final_q < -0.1);
}

#[test]
fn prescribed_temperature_retains_momentum_and_uses_the_original_shift() {
    let mut film = TransientCarreauFilm::new(law(0.4), 300.0, SPEC, 16).unwrap();
    let drive = FilmDrive {
        acceleration_m_s2: 0.5,
        ..FilmDrive::default()
    };
    film.step(
        0.2,
        drive,
        LIMITS,
        TransientCarreauControls::default(),
        || false,
    )
    .unwrap();
    let before = film.velocities_m_s().to_vec();
    let clock = film.accepted_time_s();
    let old_scale = film.zero_shear_diffusion_time_s();
    film.set_temperature(330.0, || false).unwrap();
    assert_eq!(film.velocities_m_s(), before);
    assert_eq!(film.accepted_time_s(), clock);
    assert!(film.zero_shear_diffusion_time_s() > old_scale);
    let mut restarted =
        TransientCarreauFilm::from_profile(law(0.4), 330.0, SPEC, before, clock).unwrap();
    let a = film
        .step(
            0.1,
            drive,
            LIMITS,
            TransientCarreauControls::default(),
            || false,
        )
        .unwrap();
    let b = restarted
        .step(
            0.1,
            drive,
            LIMITS,
            TransientCarreauControls::default(),
            || false,
        )
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(film.velocities_m_s(), restarted.velocities_m_s());
    assert!(film.set_temperature(400.0, || false).is_err());
    assert!(film.set_temperature(310.0, || true).is_err());
    assert_eq!(film.material().temperature_k(), 330.0);
}

#[test]
fn midpoint_refines_at_second_order_for_nonlinear_relaxation() {
    let cells = 12;
    let initial: Vec<_> = (0..cells)
        .map(|i| 0.6 * (core::f64::consts::FRAC_PI_2 * (i as f64 + 0.5) / cells as f64).sin())
        .collect();
    let advance = |steps: usize| {
        let mut film =
            TransientCarreauFilm::from_profile(law(0.4), 300.0, SPEC, initial.clone(), 0.0)
                .unwrap();
        for _ in 0..steps {
            let report = film
                .step(
                    0.2 / steps as f64,
                    FilmDrive::default(),
                    LIMITS,
                    TransientCarreauControls::default(),
                    || false,
                )
                .unwrap();
            near(
                report.flow.kinetic_energy_change_j_m2 + report.flow.viscous_heat_j_m2,
                0.0,
                1e-11,
            );
        }
        film.velocities_m_s().to_vec()
    };
    let reference = advance(1024);
    let errors: Vec<_> = [32, 64, 128]
        .into_iter()
        .map(|steps| {
            advance(steps)
                .iter()
                .zip(&reference)
                .map(|(u, v)| (u - v).powi(2))
                .sum::<f64>()
                .sqrt()
        })
        .collect();
    assert!(errors[0] / errors[1] > 3.8, "{errors:?}");
    assert!(errors[1] / errors[2] > 3.8, "{errors:?}");
}

#[test]
fn source_work_and_late_cancellation_refuse_without_changing_accepted_state() {
    let controls = TransientCarreauControls::default();
    let drive = FilmDrive {
        acceleration_m_s2: 0.8,
        surface_shear_pa: 0.1,
        ..FilmDrive::default()
    };
    let make = || TransientCarreauFilm::new(law(0.4), 300.0, SPEC, 16).unwrap();
    let mut reference = make();
    let mut polls = 0;
    let expected = reference
        .step(0.1, drive, LIMITS, controls, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut victim = make();
    let before = victim.velocities_m_s().to_vec();
    assert!(
        victim
            .step(
                0.1,
                FilmDrive {
                    surface_shear_pa: 1e6,
                    ..drive
                },
                LIMITS,
                controls,
                || false
            )
            .is_err()
    );
    assert!(matches!(
        victim.step(
            0.1,
            drive,
            LIMITS,
            TransientCarreauControls {
                maximum_material_evaluations: 1,
                ..controls
            },
            || false
        ),
        Err(TransientCarreauError::MaterialBudget | TransientCarreauError::Material(_))
    ));
    assert!(
        victim
            .step(
                0.1,
                drive,
                LIMITS,
                TransientCarreauControls {
                    maximum_newton_iterations: 1,
                    ..controls
                },
                || false
            )
            .is_err()
    );
    let mut calls = 0;
    assert!(matches!(
        victim.step(0.1, drive, LIMITS, controls, || {
            calls += 1;
            calls == polls
        }),
        Err(TransientCarreauError::Flow(TransientFilmError::Cancelled))
    ));
    assert_eq!(victim.velocities_m_s(), before);
    assert_eq!(victim.accepted_time_s(), 0.0);
    let accepted = victim.step(0.1, drive, LIMITS, controls, || false).unwrap();
    assert_eq!(accepted, expected);
    assert_eq!(victim.velocities_m_s(), reference.velocities_m_s());
}

#[test]
fn nonlinear_uniform_couette_state_and_rest_are_retained() {
    let cells = 24;
    let wall = -0.2;
    let rate = 0.7;
    let source = law(0.4);
    let traction = source
        .evaluate(300.0, rate)
        .unwrap()
        .dynamic_viscosity_pa_s()
        * rate;
    let initial: Vec<_> = (0..cells)
        .map(|i| wall + rate * (i as f64 + 0.5) / cells as f64)
        .collect();
    let mut film =
        TransientCarreauFilm::from_profile(source, 300.0, SPEC, initial.clone(), 0.0).unwrap();
    let drive = FilmDrive {
        wall_velocity_m_s: wall,
        surface_shear_pa: traction,
        ..FilmDrive::default()
    };
    let report = film
        .step(
            0.3,
            drive,
            LIMITS,
            TransientCarreauControls::default(),
            || false,
        )
        .unwrap();
    for (&a, &b) in film.velocities_m_s().iter().zip(&initial) {
        near(a, b, 2e-12);
    }
    near(report.flow.viscous_heat_j_m2, 0.3 * traction * rate, 1e-12);
    near(
        report.flow.body_work_j_m2 + report.flow.wall_work_j_m2 + report.flow.surface_work_j_m2,
        report.flow.viscous_heat_j_m2,
        2e-12,
    );
    let mut rest = TransientCarreauFilm::new(source, 300.0, SPEC, cells).unwrap();
    let zero = rest
        .step(
            1.0,
            FilmDrive::default(),
            LIMITS,
            TransientCarreauControls::default(),
            || false,
        )
        .unwrap();
    assert_eq!(zero.newton_iterations, 0);
    assert_eq!(zero.flow.viscous_heat_j_m2, 0.0);
    assert!(rest.velocities_m_s().iter().all(|u| *u == 0.0));
}

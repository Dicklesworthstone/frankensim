use super::*;

fn newtonian_time(x0: f64, x1: f64, pressure: f64, gradient: f64) -> f64 {
    let mobility = 0.01_f64.powi(2) / (8.0 * 0.1);
    (-(x1 - x0) / gradient
        - pressure / gradient.powi(2)
            * ((pressure - gradient * x1) / (pressure - gradient * x0)).ln())
        / mobility
}

#[test]
fn newtonian_rise_and_drainage_match_analytic_time_mass_and_energy() {
    let tube = model(0.1, false, 0.1, 0.0);
    for (x0, x1) in [(0.1, 0.4), (0.8, 0.6)] {
        let step = tube
            .step_inclined(
                CarreauMeniscusInput {
                    wet_length_m: x0,
                    applied_pressure_pa: 80.0,
                    duration_s: newtonian_time(x0, x1, 100.0, 200.0),
                },
                1000.0,
                0.2,
                controls(),
            )
            .unwrap();
        let m = step.meniscus;
        near(m.wet_length_m, x1, 2e-9);
        assert_eq!(m.end, CarreauMeniscusEnd::Duration);
        near(
            step.transferred_mass_kg,
            1000.0 * PI * 0.01_f64.powi(2) * (x1 - x0),
            2e-8,
        );
        near(step.hydrostatic_pressure_pa, 200.0 * x1, 2e-9);
        let gravity_energy = PI * 0.01_f64.powi(2) * 100.0 * (x1 * x1 - x0 * x0);
        near(m.gravitational_energy_change_j, gravity_energy, 2e-8);
        near(
            m.viscous_dissipation_j,
            m.reservoir_work_j - m.wall_energy_change_j - m.gravitational_energy_change_j,
            2e-9,
        );
        assert!(m.viscous_dissipation_j > 0.0);
        assert_eq!(m.transferred_volume_m3.signum(), (x1 - x0).signum());
        assert!(m.constitutive_work_defect_j.abs() < 1e-11);
    }
}

#[test]
fn equilibrium_is_approached_from_either_side_without_finite_crossing() {
    let tube = model(0.1, false, 0.0, 0.0);
    let bounded = CarreauMeniscusControls {
        maximum_profile_evaluations: 32768,
        ..controls()
    };
    for x0 in [0.1, 0.8] {
        let mut previous_gap = (x0 - 0.5_f64).abs();
        for relaxation_times in [2.0, 8.0, 20.0] {
            let duration_s = 20.0 * relaxation_times;
            let step = tube
                .step_inclined(
                    CarreauMeniscusInput {
                        wet_length_m: x0,
                        applied_pressure_pa: 100.0,
                        duration_s,
                    },
                    1000.0,
                    0.2,
                    bounded,
                )
                .unwrap()
                .meniscus;
            let gap = step.wet_length_m - 0.5;
            assert_eq!(gap.signum(), (x0 - 0.5).signum());
            assert!(gap.abs() < previous_gap);
            previous_gap = gap.abs();
            assert_eq!(step.end, CarreauMeniscusEnd::Duration);
            // The logarithmic oracle becomes condition-limited when converting
            // the final tiny gap to a rounded physical position.
            near(
                newtonian_time(x0, step.wet_length_m, 100.0, 200.0),
                duration_s,
                2e-7,
            );
            assert!(step.viscous_dissipation_j > 0.0);
            assert!(step.constitutive_work_defect_j.abs() < 1e-10);
        }
    }
    let stationary = tube
        .step_inclined(
            CarreauMeniscusInput {
                wet_length_m: 0.5,
                applied_pressure_pa: 100.0,
                duration_s: 1e12,
            },
            1000.0,
            0.2,
            controls(),
        )
        .unwrap()
        .meniscus;
    assert_eq!(stationary.end, CarreauMeniscusEnd::Stationary);
    assert_eq!(stationary.wet_length_m, 0.5);
    assert_eq!(stationary.viscous_dissipation_j, 0.0);
    assert_eq!(stationary.transferred_volume_m3, 0.0);
}

#[test]
fn pure_gravity_has_constant_newtonian_speed_and_exact_boundary_events() {
    let tube = model(0.1, false, 0.0, 0.0);
    let speed = 0.01_f64.powi(2) * 100.0 / (8.0 * 0.1);
    for gravity in [-0.1_f64, 0.1] {
        let direction = -gravity.signum();
        let step = tube
            .step_inclined(
                CarreauMeniscusInput {
                    wet_length_m: 0.2,
                    applied_pressure_pa: 0.0,
                    duration_s: 5.0,
                },
                1000.0,
                gravity,
                controls(),
            )
            .unwrap()
            .meniscus;
        near(step.wet_length_m, 0.2 + direction * speed * 5.0, 2e-10);
        near(
            step.viscous_dissipation_j,
            -step.gravitational_energy_change_j,
            2e-10,
        );
        assert_eq!(step.reservoir_work_j, 0.0);
        let x1 = if gravity < 0.0 { 1.0 } else { 0.02 };
        let arrival = (x1 - 0.2_f64).abs() / speed;
        let event = tube
            .step_inclined(
                CarreauMeniscusInput {
                    wet_length_m: 0.2,
                    applied_pressure_pa: 0.0,
                    duration_s: 2.0 * arrival,
                },
                1000.0,
                gravity,
                controls(),
            )
            .unwrap()
            .meniscus;
        assert_eq!(event.wet_length_m, x1);
        near(event.elapsed_s, arrival, 2e-12);
        assert_eq!(
            event.end,
            if gravity < 0.0 {
                CarreauMeniscusEnd::Full
            } else {
                CarreauMeniscusEnd::MinimumWetLength
            }
        );
        near(
            event.viscous_dissipation_j,
            -event.gravitational_energy_change_j,
            2e-12,
        );
    }
}

#[test]
fn nonlinear_gravity_matches_independent_radial_oracle_and_reversal() {
    let tube = model(0.4, true, 0.0, 0.0);
    let target = 0.31;
    let duration_s = PI
        * 0.01_f64.powi(2)
        * simpson(
            |x| 1.0 / independent_flow(x, 100.0 - 100.0 * x),
            0.2,
            target,
            256,
        );
    let forward = tube
        .step_inclined(
            CarreauMeniscusInput {
                wet_length_m: 0.2,
                applied_pressure_pa: 100.0,
                duration_s,
            },
            1000.0,
            0.1,
            controls(),
        )
        .unwrap()
        .meniscus;
    near(forward.wet_length_m, target, 2e-8);
    near(
        forward.endpoint_profile.flow().volumetric_flow_m3_s,
        independent_flow(target, 100.0 - 100.0 * target),
        1e-7,
    );
    near(
        forward.viscous_dissipation_j,
        forward.reservoir_work_j - forward.gravitational_energy_change_j,
        2e-9,
    );
    let reverse = tube
        .step_inclined(
            CarreauMeniscusInput {
                wet_length_m: forward.wet_length_m,
                applied_pressure_pa: -100.0,
                duration_s,
            },
            1000.0,
            -0.1,
            controls(),
        )
        .unwrap()
        .meniscus;
    near(reverse.wet_length_m, 0.2, 3e-8);
    near(
        reverse.viscous_dissipation_j,
        forward.viscous_dissipation_j,
        2e-8,
    );
    assert!(reverse.endpoint_profile.flow().volumetric_flow_m3_s < 0.0);
}

#[test]
fn gravity_steps_compose_and_horizontal_delegation_is_unchanged() {
    let tube = model(0.4, true, 0.072, 0.2);
    let input = CarreauMeniscusInput {
        wet_length_m: 0.2,
        applied_pressure_pa: 70.0,
        duration_s: 0.4,
    };
    assert_eq!(
        tube.step_horizontal(input, controls()).unwrap(),
        tube.step_inclined(input, 1000.0, 0.0, controls())
            .unwrap()
            .meniscus
    );
    let whole = tube
        .step_inclined(input, 1000.0, 0.1, controls())
        .unwrap()
        .meniscus;
    let first = tube
        .step_inclined(
            CarreauMeniscusInput {
                duration_s: 0.2,
                ..input
            },
            1000.0,
            0.1,
            controls(),
        )
        .unwrap()
        .meniscus;
    let second = tube
        .step_inclined(
            CarreauMeniscusInput {
                wet_length_m: first.wet_length_m,
                duration_s: 0.2,
                ..input
            },
            1000.0,
            0.1,
            controls(),
        )
        .unwrap()
        .meniscus;
    near(second.wet_length_m, whole.wet_length_m, 2e-9);
    near(
        first.viscous_dissipation_j + second.viscous_dissipation_j,
        whole.viscous_dissipation_j,
        2e-8,
    );
    near(
        first.gravitational_energy_change_j + second.gravitational_energy_change_j,
        whole.gravitational_energy_change_j,
        2e-8,
    );
}

#[test]
fn inclined_inputs_support_and_work_limits_are_checked_before_acceptance() {
    let tube = model(0.1, false, 0.0, 0.0);
    let input = CarreauMeniscusInput {
        wet_length_m: 0.2,
        applied_pressure_pa: 100.0,
        duration_s: 1.0,
    };
    for (density, gravity) in [
        (0.0, 0.1),
        (-1.0, 0.1),
        (f64::NAN, 0.1),
        (1000.0, f64::INFINITY),
    ] {
        assert!(matches!(
            tube.step_inclined(input, density, gravity, controls()),
            Err(CarreauMeniscusError::Input(_))
        ));
    }
    assert_eq!(
        tube.step_inclined(
            input,
            1000.0,
            0.2,
            CarreauMeniscusControls {
                maximum_profile_evaluations: 1,
                ..controls()
            }
        ),
        Err(CarreauMeniscusError::ProfileBudget)
    );
    // A weak or stationary initial load does not admit unsupported shear at
    // the declared far endpoint: the whole operating interval is checked.
    let limited = CarreauCapillaryTube::new(
        LaminarCarreauPipe::new(law(0.1, false), 300.0, 0.01, 1.0, 1.0).unwrap(),
        0.02,
        0.0,
        0.0,
    )
    .unwrap();
    for pressure in [20.1, 25.0] {
        assert!(matches!(
            limited.step_inclined(
                CarreauMeniscusInput {
                    wet_length_m: if pressure == 25.0 { 0.25 } else { 0.2 },
                    applied_pressure_pa: pressure,
                    duration_s: 0.1,
                },
                1000.0,
                0.1,
                controls()
            ),
            Err(CarreauMeniscusError::Material(_))
        ));
    }
    assert_eq!(
        tube.step_inclined(
            CarreauMeniscusInput {
                duration_s: 1e12,
                ..input
            },
            1000.0,
            0.2,
            CarreauMeniscusControls {
                maximum_profile_evaluations: 32768,
                ..controls()
            }
        ),
        Err(CarreauMeniscusError::Unrepresentable)
    );
}

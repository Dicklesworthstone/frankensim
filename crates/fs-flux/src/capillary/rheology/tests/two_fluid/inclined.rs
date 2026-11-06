use super::*;

#[test]
fn equal_viscosities_resolve_exponential_equilibrium_and_each_phase_heat() {
    let tube = pair(0.1, false, 0.1, false);
    let area = PI * 0.01_f64.powi(2);
    // Both column weights: drive = 150 - 0.2*(1000*x+500*(1-x)) = 50-100*x.
    let relaxation_s = 8.0 * 0.1 / (0.01_f64.powi(2) * 100.0);
    for (x0, multiples) in [(0.2, 12.0), (0.8, 3.0)] {
        let input = CarreauMeniscusInput {
            wet_length_m: x0,
            applied_pressure_pa: 150.0,
            duration_s: multiples * relaxation_s,
        };
        let step = tube
            .step_inclined(input, [1000.0, 500.0], 0.2, paired_controls())
            .unwrap();
        let x1 = 0.5 + (x0 - 0.5) * (-multiples).exp();
        near(step.meniscus.wet_length_m, x1, 3e-9);
        near(
            (step.meniscus.wet_length_m - 0.5).abs(),
            (x1 - 0.5).abs(),
            2e-6,
        );
        assert_eq!(step.meniscus.end, CarreauTwoFluidEnd::Duration);
        assert_eq!(
            (step.meniscus.wet_length_m - 0.5).signum(),
            (x0 - 0.5).signum()
        );
        let heat1 =
            area * (25.0 * (x1 * x1 - x0 * x0) - 100.0 / 3.0 * (x1 * x1 * x1 - x0 * x0 * x0));
        let heat_total = area * (50.0 * (x1 - x0) - 50.0 * (x1 * x1 - x0 * x0));
        near(step.meniscus.phase_viscous_dissipation_j[0], heat1, 3e-8);
        near(
            step.meniscus.phase_viscous_dissipation_j[1],
            heat_total - heat1,
            3e-8,
        );
        near(
            step.gravitational_work_j,
            step.stored_potential_energy_change_j + step.reservoir_potential_energy_out_j,
            2e-15,
        );
        near(
            step.phase_mass_changes_kg[0],
            1000.0 * area * (x1 - x0),
            3e-9,
        );
        near(
            step.phase_mass_changes_kg[1],
            -500.0 * area * (x1 - x0),
            3e-9,
        );
        assert!(step.meniscus.constitutive_work_defect_j.abs() < 3e-11);
        assert!(
            step.meniscus.profile_evaluations
                < paired_controls().meniscus.maximum_profile_evaluations
        );
    }
}

fn newtonian_time(x0: f64, x1: f64, pressure: f64, mu: [f64; 2], rho: [f64; 2], g: f64) -> f64 {
    simpson(
        |x| {
            let resistance = mu[0] * x + mu[1] * (1.0 - x);
            let drive = pressure - g * (rho[0] * x + rho[1] * (1.0 - x));
            8.0 * resistance / (0.01_f64.powi(2) * drive)
        },
        x0,
        x1,
        1024,
    )
}

#[test]
fn unequal_columns_resolve_accelerating_flow_drainage_and_arrival() {
    let tube = pair(0.1, false, 0.3, false);
    for (x0, x1, pressure, arrival) in [
        (0.2, 0.65, 200.0, false),
        (0.75, 0.3, 20.0, false),
        (0.2, 0.98, 200.0, true),
    ] {
        let duration = newtonian_time(x0, x1, pressure, [0.1, 0.3], [500.0, 1000.0], 0.1);
        let step = tube
            .step_inclined(
                CarreauMeniscusInput {
                    wet_length_m: x0,
                    applied_pressure_pa: pressure,
                    duration_s: if arrival { duration * 2.0 } else { duration },
                },
                [500.0, 1000.0],
                0.1,
                paired_controls(),
            )
            .unwrap();
        near(step.meniscus.wet_length_m, x1, 3e-9);
        near(step.meniscus.elapsed_s, duration, 3e-9);
        assert_eq!(
            step.meniscus.end,
            if arrival {
                CarreauTwoFluidEnd::MinimumDisplacedLength
            } else {
                CarreauTwoFluidEnd::Duration
            }
        );
        let expected_flow =
            PI * 0.01_f64.powi(4) * (pressure - 0.1 * (500.0 * x1 + 1000.0 * (1.0 - x1)))
                / (8.0 * (0.1 * x1 + 0.3 * (1.0 - x1)));
        near(
            step.meniscus.endpoint.volumetric_flow_m3_s,
            expected_flow,
            2e-9,
        );
        assert!(step.meniscus.constitutive_work_defect_j.abs() < 3e-11);
    }
}

#[test]
fn zero_gravity_and_equal_density_recover_horizontal_original_profiles() {
    let tube = pair(0.4, true, 0.2, false);
    let input = CarreauMeniscusInput {
        wet_length_m: 0.2,
        applied_pressure_pa: 100.0,
        duration_s: 3.0,
    };
    let horizontal = tube.step_horizontal(input, paired_controls()).unwrap();
    let zero = tube
        .step_inclined(input, [1000.0, 800.0], 0.0, paired_controls())
        .unwrap();
    assert_eq!(zero.meniscus, horizontal);
    assert_eq!(zero.gravitational_work_j, 0.0);
    let inclined = tube
        .step_inclined(input, [800.0, 800.0], 0.05, paired_controls())
        .unwrap();
    let shifted = tube
        .step_horizontal(
            CarreauMeniscusInput {
                applied_pressure_pa: 60.0,
                ..input
            },
            paired_controls(),
        )
        .unwrap();
    assert_eq!(inclined.meniscus.wet_length_m, shifted.wet_length_m);
    assert_eq!(inclined.meniscus.endpoint, shifted.endpoint);
    assert_eq!(
        inclined.meniscus.phase_viscous_dissipation_j,
        shifted.phase_viscous_dissipation_j
    );
    assert_eq!(inclined.stored_potential_energy_change_j, 0.0);
    near(
        inclined.gravitational_work_j,
        inclined.reservoir_potential_energy_out_j,
        1e-15,
    );
}

fn mixed_at(x: f64) -> (f64, f64) {
    let total = 100.0 - 0.02 * (900.0 * x + 700.0 * (1.0 - x));
    let (mut lo, mut hi) = (0.0, total);
    for _ in 0..48 {
        let p1 = 0.5 * (lo + hi);
        let q1 = independent_flow(x, p1);
        let q2 = PI * 0.01_f64.powi(4) * (total - p1) / (8.0 * 0.2 * (1.0 - x));
        if q1 > q2 {
            hi = p1;
        } else {
            lo = p1;
        }
    }
    let p1 = 0.5 * (lo + hi);
    (independent_flow(x, p1), p1)
}

#[test]
fn nonlinear_gravity_matches_independent_stress_inverse_and_pressure_partition() {
    let tube = pair(0.4, true, 0.2, false);
    let area = PI * 0.01_f64.powi(2);
    let duration = area * simpson(|x| 1.0 / mixed_at(x).0, 0.2, 0.4, 64);
    let step = tube
        .step_inclined(
            CarreauMeniscusInput {
                wet_length_m: 0.2,
                applied_pressure_pa: 100.0,
                duration_s: duration,
            },
            [900.0, 700.0],
            0.02,
            paired_controls(),
        )
        .unwrap();
    near(step.meniscus.wet_length_m, 0.4, 4e-8);
    let (q, p1) = mixed_at(0.4);
    near(step.meniscus.endpoint.volumetric_flow_m3_s, q, 1e-7);
    near(
        step.meniscus.endpoint.phase_profiles[0]
            .flow()
            .pressure_drop_pa,
        p1,
        1e-7,
    );
    near(
        step.meniscus.phase_viscous_dissipation_j[0],
        area * simpson(|x| mixed_at(x).1, 0.2, 0.4, 64),
        1e-7,
    );
    assert!(step.meniscus.constitutive_work_defect_j.abs() < 3e-10);
}

#[test]
fn reversing_axis_and_exchanging_fluids_preserves_both_physical_heats() {
    let input = CarreauMeniscusInput {
        wet_length_m: 0.2,
        applied_pressure_pa: 150.0,
        duration_s: 15.0,
    };
    let forward = pair(0.1, false, 0.3, false)
        .step_inclined(input, [1000.0, 500.0], 0.2, paired_controls())
        .unwrap();
    let reverse = pair(0.3, false, 0.1, false)
        .step_inclined(
            CarreauMeniscusInput {
                wet_length_m: 0.8,
                applied_pressure_pa: -150.0,
                ..input
            },
            [500.0, 1000.0],
            -0.2,
            paired_controls(),
        )
        .unwrap();
    near(
        reverse.meniscus.wet_length_m,
        1.0 - forward.meniscus.wet_length_m,
        3e-9,
    );
    near(
        reverse.gravitational_work_j,
        forward.gravitational_work_j,
        3e-9,
    );
    for i in 0..2 {
        near(
            reverse.meniscus.phase_viscous_dissipation_j[i],
            forward.meniscus.phase_viscous_dissipation_j[1 - i],
            3e-8,
        );
    }
}

#[test]
fn stationary_support_envelope_and_unresolved_equilibrium_are_explicit() {
    let tube = pair(0.1, false, 0.1, false);
    let equilibrium = CarreauMeniscusInput {
        wet_length_m: 0.5,
        applied_pressure_pa: 150.0,
        duration_s: 100.0,
    };
    let stationary = tube
        .step_inclined(equilibrium, [1000.0, 500.0], 0.2, paired_controls())
        .unwrap();
    assert_eq!(stationary.meniscus.end, CarreauTwoFluidEnd::Stationary);
    assert_eq!(stationary.meniscus.viscous_dissipation_j, 0.0);
    assert_eq!(stationary.gravitational_work_j, 0.0);
    assert_eq!(
        tube.step_inclined(
            CarreauMeniscusInput {
                wet_length_m: 0.2,
                duration_s: 8000.0,
                ..equilibrium
            },
            [1000.0, 500.0],
            0.2,
            paired_controls()
        ),
        Err(CarreauMeniscusError::Unrepresentable)
    );
    assert!(
        tube.step_inclined(equilibrium, [0.0, 500.0], 0.2, paired_controls())
            .is_err()
    );
    let limited = LaminarCarreauPipe::new(law(0.1, false), 300.0, 0.01, 1.0, 1.0).unwrap();
    let limited = CarreauTwoFluidCapillary::new(
        CarreauCapillaryTube::new(limited, 0.02, 0.0, 0.0).unwrap(),
        limited,
        0.02,
    )
    .unwrap();
    // Initial drive is 10.4 Pa, inside the 20 Pa source limit. The other
    // declared endpoint requires 29.6 Pa, so even a zero-duration step refuses.
    assert!(matches!(
        limited.step_inclined(
            CarreauMeniscusInput {
                wet_length_m: 0.02,
                applied_pressure_pa: 40.0,
                duration_s: 0.0,
            },
            [1.0, 3.0],
            10.0,
            paired_controls()
        ),
        Err(CarreauMeniscusError::Material(_))
    ));
}

#[test]
fn inclined_late_cancellation_retains_identical_retry() {
    let tube = pair(0.1, false, 0.3, false);
    let input = CarreauMeniscusInput {
        wet_length_m: 0.2,
        applied_pressure_pa: 150.0,
        duration_s: 15.0,
    };
    let mut calls = 0;
    let expected = tube
        .step_inclined_with_check(input, [1000.0, 500.0], 0.2, paired_controls(), || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    let mut attempted = 0;
    let error = tube
        .step_inclined_with_check(input, [1000.0, 500.0], 0.2, paired_controls(), || {
            attempted += 1;
            if attempted == calls {
                Err(CarreauMeniscusError::Input("cancelled"))
            } else {
                Ok(())
            }
        })
        .unwrap_err();
    assert_eq!(error, CarreauMeniscusError::Input("cancelled"));
    assert_eq!(attempted, calls);
    assert_eq!(
        tube.step_inclined(input, [1000.0, 500.0], 0.2, paired_controls())
            .unwrap(),
        expected
    );
}

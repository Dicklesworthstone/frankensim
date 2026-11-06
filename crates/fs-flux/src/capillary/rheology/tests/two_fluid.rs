use super::super::two_fluid::{
    CarreauTwoFluidCapillary, CarreauTwoFluidControls, CarreauTwoFluidEnd,
};
use super::*;

fn pair(mu1: f64, thin1: bool, mu2: f64, thin2: bool) -> CarreauTwoFluidCapillary {
    CarreauTwoFluidCapillary::new(
        model(mu1, thin1, 0.0, 0.0),
        LaminarCarreauPipe::new(law(mu2, thin2), 300.0, 0.01, 1.0, 1e6).unwrap(),
        0.02,
    )
    .unwrap()
}

fn paired_controls() -> CarreauTwoFluidControls {
    CarreauTwoFluidControls {
        meniscus: controls(),
        relative_pressure_tolerance: 1e-10,
        ..CarreauTwoFluidControls::default()
    }
}

fn exact_time(x0: f64, x1: f64, pressure: f64, mu1: f64, mu2: f64) -> f64 {
    8.0 / (0.01_f64.powi(2) * pressure)
        * (mu2 * (x1 - x0) + 0.5 * (mu1 - mu2) * (x1 * x1 - x0 * x0))
}

#[test]
fn two_newtonian_phases_match_independent_travel_and_separate_heat() {
    let tube = pair(0.1, false, 0.3, false);
    for (x0, x1, pressure) in [(0.2, 0.65, 100.0), (0.75, 0.3, -100.0)] {
        let step = tube
            .step_horizontal(
                CarreauMeniscusInput {
                    wet_length_m: x0,
                    applied_pressure_pa: pressure,
                    duration_s: exact_time(x0, x1, pressure, 0.1, 0.3),
                },
                paired_controls(),
            )
            .unwrap();
        near(step.wet_length_m, x1, 2e-9);
        assert_eq!(step.end, CarreauTwoFluidEnd::Duration);
        let resistance = 0.1 * x1 + 0.3 * (1.0 - x1);
        let q = PI * 0.01_f64.powi(4) * pressure / (8.0 * resistance);
        near(step.endpoint.volumetric_flow_m3_s, q, 2e-9);
        let area = PI * 0.01_f64.powi(2);
        let heat1 =
            area * pressure * simpson(|x| 0.1 * x / (0.1 * x + 0.3 * (1.0 - x)), x0, x1, 512);
        let heat2 = area * pressure * (x1 - x0) - heat1;
        near(step.phase_viscous_dissipation_j[0], heat1, 2e-8);
        near(step.phase_viscous_dissipation_j[1], heat2, 2e-8);
        near(step.viscous_dissipation_j, step.reservoir_work_j, 2e-10);
        assert!(step.constitutive_work_defect_j.abs() < 1e-11);
        assert_eq!(
            step.endpoint.phase_profiles[0].pipe().length_m(),
            step.wet_length_m
        );
        near(
            step.endpoint.phase_profiles[1].pipe().length_m(),
            1.0 - x1,
            2e-9,
        );
        for residual in step.endpoint.phase_flow_residuals_m3_s {
            assert!(residual.abs() <= 2e-12 * q.abs());
        }
    }
}

#[test]
fn identical_nonlinear_fluids_have_constant_interface_speed_and_resolved_heat() {
    let tube = pair(0.4, true, 0.4, true);
    let q = independent_flow(1.0, 100.0);
    let duration_s = PI * 0.01_f64.powi(2) * 0.3 / q;
    let step = tube
        .step_horizontal(
            CarreauMeniscusInput {
                wet_length_m: 0.2,
                applied_pressure_pa: 100.0,
                duration_s,
            },
            paired_controls(),
        )
        .unwrap();
    near(step.wet_length_m, 0.5, 2e-8);
    near(step.endpoint.volumetric_flow_m3_s, q, 1e-7);
    let heat1 = PI * 0.01_f64.powi(2) * 100.0 * (0.5_f64.powi(2) - 0.2_f64.powi(2)) * 0.5;
    near(step.phase_viscous_dissipation_j[0], heat1, 2e-8);
    near(
        step.phase_viscous_dissipation_j[1],
        step.reservoir_work_j - heat1,
        2e-8,
    );
    let reverse = tube
        .step_horizontal(
            CarreauMeniscusInput {
                wet_length_m: step.wet_length_m,
                applied_pressure_pa: -100.0,
                duration_s,
            },
            paired_controls(),
        )
        .unwrap();
    near(reverse.wet_length_m, 0.2, 3e-8);
    for i in 0..2 {
        near(
            reverse.phase_viscous_dissipation_j[i],
            step.phase_viscous_dissipation_j[i],
            3e-8,
        );
    }
}

fn independent_mixed(x: f64) -> (f64, f64) {
    let (mut low, mut high) = (0.0, 100.0);
    for _ in 0..48 {
        let p1 = 0.5 * (low + high);
        let q1 = independent_flow(x, p1);
        let q2 = PI * 0.01_f64.powi(4) * (100.0 - p1) / (8.0 * 0.2 * (1.0 - x));
        if q1 > q2 {
            high = p1;
        } else {
            low = p1;
        }
    }
    let p1 = 0.5 * (low + high);
    (independent_flow(x, p1), p1)
}

#[test]
fn unequal_laws_match_independent_pressure_partition_and_travel() {
    let tube = pair(0.4, true, 0.2, false);
    let area = PI * 0.01_f64.powi(2);
    let duration_s = area * simpson(|x| 1.0 / independent_mixed(x).0, 0.2, 0.4, 64);
    let step = tube
        .step_horizontal(
            CarreauMeniscusInput {
                wet_length_m: 0.2,
                applied_pressure_pa: 100.0,
                duration_s,
            },
            paired_controls(),
        )
        .unwrap();
    near(step.wet_length_m, 0.4, 3e-8);
    let (q, p1) = independent_mixed(0.4);
    near(step.endpoint.volumetric_flow_m3_s, q, 1e-7);
    near(
        step.endpoint.phase_profiles[0].flow().pressure_drop_pa,
        p1,
        1e-7,
    );
    near(
        step.endpoint.phase_profiles[1].flow().pressure_drop_pa,
        100.0 - p1,
        1e-7,
    );
    let heat1 = area * simpson(|x| independent_mixed(x).1, 0.2, 0.4, 64);
    near(step.phase_viscous_dissipation_j[0], heat1, 1e-7);
    near(
        step.phase_viscous_dissipation_j[1],
        area * 100.0 * 0.2 - heat1,
        1e-7,
    );
    assert!(step.constitutive_work_defect_j.abs() < 2e-10);
}

#[test]
fn phase_exchange_and_temporal_subdivision_preserve_work() {
    let tube = pair(0.1, false, 0.3, false);
    let swapped = pair(0.3, false, 0.1, false);
    let input = CarreauMeniscusInput {
        wet_length_m: 0.2,
        applied_pressure_pa: 100.0,
        duration_s: 10.0,
    };
    let whole = tube.step_horizontal(input, paired_controls()).unwrap();
    let mirrored = swapped
        .step_horizontal(
            CarreauMeniscusInput {
                wet_length_m: 0.8,
                applied_pressure_pa: -100.0,
                duration_s: input.duration_s,
            },
            paired_controls(),
        )
        .unwrap();
    near(mirrored.wet_length_m, 1.0 - whole.wet_length_m, 2e-9);
    near(
        mirrored.phase_viscous_dissipation_j[0],
        whole.phase_viscous_dissipation_j[1],
        2e-8,
    );
    near(
        mirrored.phase_viscous_dissipation_j[1],
        whole.phase_viscous_dissipation_j[0],
        2e-8,
    );
    let first = tube
        .step_horizontal(
            CarreauMeniscusInput {
                duration_s: 5.0,
                ..input
            },
            paired_controls(),
        )
        .unwrap();
    let second = tube
        .step_horizontal(
            CarreauMeniscusInput {
                wet_length_m: first.wet_length_m,
                duration_s: 5.0,
                ..input
            },
            paired_controls(),
        )
        .unwrap();
    near(second.wet_length_m, whole.wet_length_m, 2e-9);
    for i in 0..2 {
        near(
            first.phase_viscous_dissipation_j[i] + second.phase_viscous_dissipation_j[i],
            whole.phase_viscous_dissipation_j[i],
            2e-8,
        );
    }
}

#[test]
fn declared_column_limits_return_arrival_time_and_stationary_balance() {
    let tube = pair(0.1, false, 0.3, false);
    for (x1, pressure, end) in [
        (0.98, 100.0, CarreauTwoFluidEnd::MinimumDisplacedLength),
        (0.02, -100.0, CarreauTwoFluidEnd::MinimumInvadingLength),
    ] {
        let arrival = exact_time(0.2, x1, pressure, 0.1, 0.3);
        let step = tube
            .step_horizontal(
                CarreauMeniscusInput {
                    wet_length_m: 0.2,
                    applied_pressure_pa: pressure,
                    duration_s: 2.0 * arrival,
                },
                paired_controls(),
            )
            .unwrap();
        assert_eq!(step.wet_length_m, x1);
        assert_eq!(step.end, end);
        near(step.elapsed_s, arrival, 2e-12);
        assert!(step.phase_viscous_dissipation_j.iter().all(|h| *h > 0.0));
    }
    let stationary = tube
        .step_horizontal(
            CarreauMeniscusInput {
                wet_length_m: 0.2,
                applied_pressure_pa: 0.0,
                duration_s: 100.0,
            },
            paired_controls(),
        )
        .unwrap();
    assert_eq!(stationary.end, CarreauTwoFluidEnd::Stationary);
    assert_eq!(stationary.transferred_volume_m3, 0.0);
    assert_eq!(stationary.phase_viscous_dissipation_j, [0.0; 2]);
}

#[test]
fn original_thermal_states_and_capillary_work_are_retained() {
    let invading = model(0.1, false, 0.072, 0.0);
    let displaced = LaminarCarreauPipe::new(law(0.3, false), 320.0, 0.01, 1.0, 1e6).unwrap();
    let tube = CarreauTwoFluidCapillary::new(invading, displaced, 0.02).unwrap();
    let step = tube
        .step_horizontal(
            CarreauMeniscusInput {
                wet_length_m: 0.2,
                applied_pressure_pa: 0.0,
                duration_s: 5.0,
            },
            CarreauTwoFluidControls::default(),
        )
        .unwrap();
    assert!(step.wet_length_m > 0.2);
    assert_eq!(step.reservoir_work_j, 0.0);
    near(step.viscous_dissipation_j, -step.wall_energy_change_j, 2e-8);
    assert_eq!(
        step.endpoint.phase_profiles[0].pipe().temperature_k(),
        300.0
    );
    assert_eq!(
        step.endpoint.phase_profiles[1].pipe().temperature_k(),
        320.0
    );
    assert_eq!(
        *step.endpoint.phase_profiles[0].pipe().law(),
        *tube.invading().pipe().law()
    );
    assert_eq!(
        *step.endpoint.phase_profiles[1].pipe().law(),
        *tube.displaced().law()
    );
}

#[test]
fn whole_interval_support_and_each_work_budget_refuse_without_mutation() {
    let tube = pair(0.4, true, 0.2, false);
    let original = tube;
    let input = CarreauMeniscusInput {
        wet_length_m: 0.2,
        applied_pressure_pa: 100.0,
        duration_s: 1.0,
    };
    assert_eq!(
        tube.step_horizontal(
            input,
            CarreauTwoFluidControls {
                meniscus: CarreauMeniscusControls {
                    maximum_profile_evaluations: 1,
                    ..controls()
                },
                ..paired_controls()
            }
        ),
        Err(CarreauMeniscusError::ProfileBudget)
    );
    assert_eq!(
        tube.step_horizontal(
            input,
            CarreauTwoFluidControls {
                maximum_pressure_iterations: 1,
                ..paired_controls()
            }
        ),
        Err(CarreauMeniscusError::SeriesPressureBudget)
    );
    let limited1 = CarreauCapillaryTube::new(
        LaminarCarreauPipe::new(law(0.1, false), 300.0, 0.01, 1.0, 1.0).unwrap(),
        0.02,
        0.0,
        0.0,
    )
    .unwrap();
    let limited = CarreauTwoFluidCapillary::new(
        limited1,
        LaminarCarreauPipe::new(law(0.3, false), 300.0, 0.01, 1.0, 1.0).unwrap(),
        0.02,
    )
    .unwrap();
    assert!(matches!(
        limited.step_horizontal(
            CarreauMeniscusInput {
                applied_pressure_pa: 40.0,
                ..input
            },
            paired_controls()
        ),
        Err(CarreauMeniscusError::Material(_))
    ));
    assert!(matches!(
        CarreauTwoFluidCapillary::new(*tube.invading(), *tube.displaced(), 0.0),
        Err(CarreauMeniscusError::Input(_))
    ));
    assert!(matches!(
        tube.step_horizontal(
            CarreauMeniscusInput {
                wet_length_m: 1.0,
                ..input
            },
            paired_controls()
        ),
        Err(CarreauMeniscusError::Input(_))
    ));
    assert_eq!(tube, original);
    assert!(
        tube.step_horizontal(input, CarreauTwoFluidControls::default())
            .is_ok()
    );
}

#[test]
fn cancellation_at_late_material_query_and_publication_preserves_exact_retry() {
    let tube = pair(0.1, false, 0.3, false);
    let input = CarreauMeniscusInput {
        wet_length_m: 0.2,
        applied_pressure_pa: 100.0,
        duration_s: 10.0,
    };
    let mut calls = 0;
    let reference = tube
        .step_horizontal_with_check(input, paired_controls(), || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    assert_eq!(calls, reference.profile_evaluations + 2);
    let cancelled = CarreauMeniscusError::Input("caller cancelled");
    for stop in [calls - 1, calls] {
        let mut seen = 0;
        let result = tube.step_horizontal_with_check(input, paired_controls(), || {
            seen += 1;
            if seen == stop {
                Err(cancelled.clone())
            } else {
                Ok(())
            }
        });
        assert_eq!(result, Err(cancelled.clone()));
        assert_eq!(
            tube.step_horizontal(input, paired_controls()).unwrap(),
            reference
        );
    }
}

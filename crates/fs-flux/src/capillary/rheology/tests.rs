use super::*;
mod inclined;
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-20),
        "{actual:.16e} != {expected:.16e}"
    );
}

fn law(mu: f64, thinning: bool) -> ThermalCarreauYasuda {
    let curve = CarreauYasuda::new(
        mu,
        if thinning { 0.0 } else { mu },
        0.8,
        0.5,
        1.0,
        [0.0, 1e6],
    )
    .unwrap();
    ThermalCarreauYasuda::new(
        curve,
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 4000.0,
            },
            300.0,
            [280.0, 350.0],
        )
        .unwrap(),
        ViscosityShiftMode::Vertical,
    )
}

fn model(mu: f64, thinning: bool, tension: f64, angle: f64) -> CarreauCapillaryTube {
    let pipe = LaminarCarreauPipe::new(law(mu, thinning), 300.0, 0.01, 1.0, 1e6).unwrap();
    CarreauCapillaryTube::new(pipe, 0.02, tension, angle).unwrap()
}

fn controls() -> CarreauMeniscusControls {
    CarreauMeniscusControls {
        pipe: CarreauPipeControls {
            relative_constitutive_tolerance: 1e-12,
            relative_quadrature_tolerance: 1e-10,
            ..CarreauPipeControls::default()
        },
        relative_time_tolerance: 1e-9,
        ..CarreauMeniscusControls::default()
    }
}

#[test]
fn newtonian_pressure_filling_matches_independent_washburn_and_work() {
    let tube = model(0.1, false, 0.0, 0.0);
    let input = CarreauMeniscusInput {
        wet_length_m: 0.2,
        applied_pressure_pa: 100.0,
        duration_s: 0.4,
    };
    let step = tube.step_horizontal(input, controls()).unwrap();
    let expected = (input.wet_length_m.powi(2)
        + 0.01_f64.powi(2) * 100.0 * input.duration_s / (4.0 * 0.1))
        .sqrt();
    near(step.wet_length_m, expected, 2e-10);
    assert_eq!(step.end, CarreauMeniscusEnd::Duration);
    assert_eq!(step.elapsed_s, input.duration_s);
    near(
        step.transferred_volume_m3,
        PI * 0.01_f64.powi(2) * (expected - 0.2),
        2e-9,
    );
    near(step.viscous_dissipation_j, step.reservoir_work_j, 2e-12);
    assert_eq!(step.wall_energy_change_j, 0.0);
    assert!(step.constitutive_work_defect_j.abs() < 1e-14);
    assert!(
        step.travel_time_residual_s.abs() <= controls().relative_time_tolerance * input.duration_s
    );
    assert_eq!(*step.endpoint_profile.pipe().law(), *tube.pipe().law());
    assert_eq!(step.endpoint_profile.pipe().length_m(), step.wet_length_m);
    assert_eq!(step.endpoint_profile.pipe().temperature_k(), 300.0);
}

#[test]
fn surface_energy_drives_filling_and_signed_pressure_drives_drainage() {
    let tube = model(0.1, false, 0.072, PI / 3.0);
    let input = CarreauMeniscusInput {
        wet_length_m: 0.2,
        applied_pressure_pa: 0.0,
        duration_s: 2.0,
    };
    let step = tube.step_horizontal(input, controls()).unwrap();
    let expected = (0.2_f64.powi(2) + 0.01 * 0.072 * 0.5 / (2.0 * 0.1) * 2.0).sqrt();
    near(step.wet_length_m, expected, 2e-10);
    assert!(step.wall_energy_change_j < 0.0);
    assert_eq!(step.reservoir_work_j, 0.0);
    near(
        step.viscous_dissipation_j,
        -step.wall_energy_change_j,
        2e-12,
    );

    let reverse = tube
        .step_horizontal(
            CarreauMeniscusInput {
                wet_length_m: step.wet_length_m,
                applied_pressure_pa: -2.0 * tube.capillary_pressure_pa(),
                duration_s: input.duration_s,
            },
            controls(),
        )
        .unwrap();
    near(reverse.wet_length_m, input.wet_length_m, 3e-10);
    assert!(reverse.transferred_volume_m3 < 0.0);
    assert!(reverse.wall_energy_change_j > 0.0);
    assert!(reverse.reservoir_work_j > reverse.viscous_dissipation_j);
    near(
        reverse.reservoir_work_j - reverse.wall_energy_change_j,
        reverse.viscous_dissipation_j,
        2e-12,
    );
}

#[test]
fn geometric_events_consume_only_arrival_time_without_dry_singularity() {
    let tube = model(0.1, false, 0.0, 0.0);
    let speed_coefficient = 0.01_f64.powi(2) * 100.0 / (8.0 * 0.1);
    for (x0, pressure, end, x1) in [
        (0.9_f64, 100.0, CarreauMeniscusEnd::Full, 1.0_f64),
        (0.2, -100.0, CarreauMeniscusEnd::MinimumWetLength, 0.02),
    ] {
        let arrival = (x1 * x1 - x0 * x0).abs() / (2.0 * speed_coefficient);
        let step = tube
            .step_horizontal(
                CarreauMeniscusInput {
                    wet_length_m: x0,
                    applied_pressure_pa: pressure,
                    duration_s: 2.0 * arrival,
                },
                controls(),
            )
            .unwrap();
        assert_eq!(step.end, end);
        assert_eq!(step.wet_length_m, x1);
        near(step.elapsed_s, arrival, 2e-12);
        assert!(step.elapsed_s < step.input.duration_s);
        assert!(step.viscous_dissipation_j > 0.0);
        let stopped = tube
            .step_horizontal(
                CarreauMeniscusInput {
                    wet_length_m: x1,
                    applied_pressure_pa: pressure,
                    duration_s: arrival,
                },
                controls(),
            )
            .unwrap();
        assert_eq!(stopped.elapsed_s, 0.0);
        assert_eq!(stopped.end, end);
        assert_eq!(stopped.transferred_volume_m3, 0.0);
    }
}

fn simpson(f: impl Fn(f64) -> f64, a: f64, b: f64, panels: u32) -> f64 {
    let dx = (b - a) / panels as f64;
    let mut sum = f(a) + f(b);
    for i in 1..panels {
        sum += if i % 2 == 0 { 2.0 } else { 4.0 } * f(a + dx * i as f64);
    }
    sum * dx / 3.0
}

// Independent algebraic inverse for mu=mu0/sqrt(1+lambda*gamma):
// s=lambda*tau/mu0, t=lambda*gamma=(s^2+s*sqrt(s^2+4))/2.
// Radial momentum then gives Q=pi R^3 integral(z^2 gamma(tau_w*z) dz).
fn independent_flow(length: f64, pressure: f64) -> f64 {
    let r = 0.01;
    let tau_w = pressure * r / (2.0 * length);
    PI * r
        * r
        * r
        * simpson(
            |z| {
                let s = 0.8 * tau_w * z / 0.4;
                let rate = 0.5 * (s * s + s * (s * s + 4.0).sqrt()) / 0.8;
                z * z * rate
            },
            0.0,
            1.0,
            512,
        )
}

#[test]
fn nonlinear_filling_matches_independent_radial_and_travel_integrals() {
    let tube = model(0.4, true, 0.0, 0.0);
    let target = 0.31;
    let elapsed =
        PI * 0.01_f64.powi(2) * simpson(|x| 1.0 / independent_flow(x, 100.0), 0.2, target, 256);
    let step = tube
        .step_horizontal(
            CarreauMeniscusInput {
                wet_length_m: 0.2,
                applied_pressure_pa: 100.0,
                duration_s: elapsed,
            },
            controls(),
        )
        .unwrap();
    near(step.wet_length_m, target, 2e-8);
    near(
        step.endpoint_profile.flow().volumetric_flow_m3_s,
        independent_flow(target, 100.0),
        1e-7,
    );
    let frozen_zero_shear =
        (0.2_f64.powi(2) + 0.01_f64.powi(2) * 100.0 * elapsed / (4.0 * 0.4)).sqrt();
    assert!(step.wet_length_m > frozen_zero_shear * 1.2);
    near(step.viscous_dissipation_j, step.reservoir_work_j, 2e-10);
    assert!(step.profile_evaluations <= controls().maximum_profile_evaluations);
    let half = tube
        .step_horizontal(
            CarreauMeniscusInput {
                duration_s: elapsed * 0.5,
                ..step.input
            },
            controls(),
        )
        .unwrap();
    let second = tube
        .step_horizontal(
            CarreauMeniscusInput {
                wet_length_m: half.wet_length_m,
                duration_s: elapsed * 0.5,
                ..step.input
            },
            controls(),
        )
        .unwrap();
    near(second.wet_length_m, step.wet_length_m, 3e-9);
    near(
        half.viscous_dissipation_j + second.viscous_dissipation_j,
        step.viscous_dissipation_j,
        1e-8,
    );
}

#[test]
fn nonlinear_reversal_and_stationary_state_retain_original_law() {
    let tube = model(0.4, true, 0.0, 0.0);
    let input = CarreauMeniscusInput {
        wet_length_m: 0.2,
        applied_pressure_pa: 80.0,
        duration_s: 0.1,
    };
    let forward = tube.step_horizontal(input, controls()).unwrap();
    let reverse = tube
        .step_horizontal(
            CarreauMeniscusInput {
                wet_length_m: forward.wet_length_m,
                applied_pressure_pa: -80.0,
                duration_s: input.duration_s,
            },
            controls(),
        )
        .unwrap();
    near(reverse.wet_length_m, input.wet_length_m, 5e-9);
    near(
        reverse.viscous_dissipation_j,
        forward.viscous_dissipation_j,
        5e-8,
    );
    assert!(reverse.endpoint_profile.flow().volumetric_flow_m3_s < 0.0);
    let resting = tube
        .step_horizontal(
            CarreauMeniscusInput {
                applied_pressure_pa: 0.0,
                duration_s: 10.0,
                ..input
            },
            controls(),
        )
        .unwrap();
    assert_eq!(resting.wet_length_m, input.wet_length_m);
    assert_eq!(resting.elapsed_s, 10.0);
    assert_eq!(resting.end, CarreauMeniscusEnd::Stationary);
    assert_eq!(resting.viscous_dissipation_j, 0.0);
    assert_eq!(resting.endpoint_profile.flow().volumetric_flow_m3_s, 0.0);
}

#[test]
fn support_and_work_budgets_refuse_without_mutating_model_or_input() {
    let tube = model(0.4, true, 0.0, 0.0);
    let saved = tube;
    let input = CarreauMeniscusInput {
        wet_length_m: 0.2,
        applied_pressure_pa: 80.0,
        duration_s: 0.1,
    };
    assert!(matches!(
        tube.step_horizontal(
            input,
            CarreauMeniscusControls {
                maximum_profile_evaluations: 1,
                ..controls()
            }
        ),
        Err(CarreauMeniscusError::ProfileBudget)
    ));
    assert!(matches!(
        tube.step_horizontal(
            input,
            CarreauMeniscusControls {
                maximum_iterations: 0,
                ..controls()
            }
        ),
        Err(CarreauMeniscusError::Input(_))
    ));
    assert!(matches!(
        tube.step_horizontal(
            CarreauMeniscusInput {
                wet_length_m: 0.0,
                ..input
            },
            controls()
        ),
        Err(CarreauMeniscusError::Input(_))
    ));
    assert!(matches!(
        tube.step_horizontal(
            CarreauMeniscusInput {
                duration_s: -1.0,
                ..input
            },
            controls()
        ),
        Err(CarreauMeniscusError::Input(_))
    ));
    let limited = LaminarCarreauPipe::new(law(0.4, true), 300.0, 0.01, 1.0, 1.0).unwrap();
    let limited = CarreauCapillaryTube::new(limited, 0.02, 0.0, 0.0).unwrap();
    assert!(matches!(
        limited.step_horizontal(input, controls()),
        Err(CarreauMeniscusError::Material(_))
    ));
    assert_eq!(tube, saved);
    assert_eq!(input.wet_length_m, 0.2);
    assert!(tube.pipe().with_length(0.0).is_err());
    assert!(tube.pipe().with_length(f64::INFINITY).is_err());
}

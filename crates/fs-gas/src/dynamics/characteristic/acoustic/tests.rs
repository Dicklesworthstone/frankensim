use crate::dynamics::tests::gas;
use crate::dynamics::{Boundaries, Boundary, Error, Euler1d, FaceFlux, StepControl};
use fs_qty::{Length, Time};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-9 * (1.0 + expected.abs()),
        "actual={actual}, expected={expected}"
    );
}

fn close_flux(actual: FaceFlux, expected: FaceFlux) {
    close(actual.mass, expected.mass);
    close(actual.momentum, expected.momentum);
    close(actual.energy, expected.energy);
}

// Independent Simpson quadrature of Cp(T)/sqrt(gamma(T)*R*T), R=1.
fn acoustic_integral(temperature: f64, linear_cp: f64) -> f64 {
    let n = 512;
    let width = (temperature - 1.0) / f64::from(n);
    let mut sum = 0.0;
    for i in 0..=n {
        let t = 1.0 + f64::from(i) * width;
        let cp = 3.5 + linear_cp * t;
        let c = ((cp / (cp - 1.0)) * t).sqrt();
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        sum += weight * cp / c;
    }
    width * sum / 3.0
}

#[test]
fn characteristic_outlet_passes_outgoing_simple_waves_with_original_cp_and_energy_datum() {
    for linear_cp in [0.0, 0.4] {
        let source = gas(linear_cp, 0.0);
        let shifted = gas(linear_cp, 1e18);
        for normal in [-1.0, 1.0] {
            for temperature in [0.95_f64, 1.0, 1.05] {
                // s-s0=3.5 ln(T)+a(T-1)-ln(p), independently of the solver.
                let pressure = temperature.powf(3.5) * (linear_cp * (temperature - 1.0)).exp();
                let density = pressure / temperature;
                let velocity = normal * (0.4 + acoustic_integral(temperature, linear_cp));
                let state = source.conserved(density, velocity, pressure).unwrap();
                let inside = source.primitive(state, || false).unwrap();
                let expected = source.face_flux(inside, 0.0).unwrap();
                let actual = source
                    .characteristic_outlet_flux(state, normal, 1.0, normal * 0.4, || false)
                    .unwrap();
                close_flux(actual, expected);
                assert_eq!(
                    actual,
                    shifted
                        .characteristic_outlet_flux(state, normal, 1.0, normal * 0.4, || false)
                        .unwrap()
                );
            }
        }
    }
}

#[test]
fn characteristic_outlet_passes_entropy_contacts_without_turning_temperature_into_sound() {
    let source = gas(0.4, 0.0);
    for normal in [-1.0, 1.0] {
        for density in [0.3, 1.0, 3.0] {
            let state = source.conserved(density, normal * 0.2, 1.0).unwrap();
            let inside = source.primitive(state, || false).unwrap();
            let actual = source
                .characteristic_outlet_flux(state, normal, 1.0, normal * 0.2, || false)
                .unwrap();
            close_flux(actual, source.face_flux(inside, 0.0).unwrap());
        }
    }
}

#[test]
fn characteristic_outlet_prescribes_only_the_incoming_gamma_characteristic() {
    let source = gas(0.0, 0.0);
    let initial_sound_speed = 1.4_f64.sqrt();
    for normal in [-1.0, 1.0] {
        for excess_velocity in [2e-8, 0.4] {
            let state = source
                .conserved(1.0, normal * (0.2 + excess_velocity), 1.0)
                .unwrap();
            // J+ = u + 5c and J- = u - 5c: only half the velocity
            // perturbation remains; the other half becomes pressure.
            let velocity = normal * (0.2 + 0.5 * excess_velocity);
            let sound_speed = initial_sound_speed + excess_velocity / 10.0;
            let temperature = (sound_speed / initial_sound_speed).powi(2);
            let expected = source
                .from_rho_u_p(temperature.powf(2.5), velocity, temperature.powf(3.5))
                .unwrap();
            close_flux(
                source
                    .characteristic_outlet_flux(state, normal, 1.0, normal * 0.2, || false)
                    .unwrap(),
                source.face_flux(expected, 0.0).unwrap(),
            );
        }
    }
}

#[test]
fn characteristic_outlet_respects_sonic_discharge_and_supersonic_causality() {
    let source = gas(0.0, 0.0);
    for normal in [-1.0, 1.0] {
        let state = source.conserved(1.0, normal * 3.0, 1.0).unwrap();
        let inside = source.primitive(state, || false).unwrap();
        for pressure in [1e-100, 1e100] {
            close_flux(
                source
                    .characteristic_outlet_flux(state, normal, pressure, normal * 0.2, || false)
                    .unwrap(),
                source.face_flux(inside, 0.0).unwrap(),
            );
        }
        let state = source.conserved(1.0, normal * 0.2, 1.0).unwrap();
        let sound_speed = (0.2 + 5.0 * 1.4_f64.sqrt()) / 6.0;
        let temperature = sound_speed * sound_speed / 1.4;
        let sonic = source
            .from_rho_u_p(
                temperature.powf(2.5),
                normal * sound_speed,
                temperature.powf(3.5),
            )
            .unwrap();
        close_flux(
            source
                .characteristic_outlet_flux(state, normal, 0.001, normal * 0.2, || false)
                .unwrap(),
            source.face_flux(sonic, 0.0).unwrap(),
        );
    }
}

#[test]
fn characteristic_outlet_refuses_backflow_invalid_inputs_and_unsupported_isentropes() {
    let source = gas(0.2, 0.0);
    let state = source.conserved(1.0, 0.2, 1.0).unwrap();
    for (normal, pressure, velocity) in [
        (0.0, 1.0, 0.2),
        (1.0, 0.0, 0.2),
        (1.0, f64::NAN, 0.2),
        (1.0, 1.0, f64::INFINITY),
        (1.0, 1.0, -0.2),
        (-1.0, 1.0, -0.2),
        (1.0, 100.0, 0.2),
    ] {
        assert!(matches!(
            source.characteristic_outlet_flux(state, normal, pressure, velocity, || false),
            Err(Error::Invalid(_))
        ));
    }
    assert!(matches!(
        source.characteristic_outlet_flux(state, 1.0, 1e-40, 0.2, || false),
        Err(Error::CaloricPath(_))
    ));
    assert!(matches!(
        source.characteristic_outlet_flux(state, 1.0, 0.9, 0.2, || true),
        Err(Error::Cancelled)
    ));
}

fn outgoing_pulse_remaining_pressure(right: Boundary) -> f64 {
    let source = gas(0.0, 0.0);
    let cells = 96;
    let states = (0..cells)
        .map(|i| {
            let x = (f64::from(i) + 0.5) / f64::from(cells);
            let pressure = 1.0 + 1e-3 * (-((x - 0.65) / 0.055).powi(2)).exp();
            let density = pressure.powf(1.0 / 1.4);
            let c = (1.4 * pressure / density).sqrt();
            let velocity = 0.3 + 5.0 * (c - 1.4_f64.sqrt());
            source.conserved(density, velocity, pressure).unwrap()
        })
        .collect();
    let mut flow = Euler1d::new(
        source,
        Length::new(1.0 / f64::from(cells)),
        1.0,
        states,
        Boundaries::Ends {
            left: Boundary::Transmissive,
            right,
        },
        || false,
    )
    .unwrap();
    while flow.time().value() < 0.6 {
        let report = flow
            .step(
                Time::new(0.6 - flow.time().value()),
                StepControl::default(),
                || false,
            )
            .unwrap();
        assert!(report.conservation_residual.mass.abs() < 1e-12);
        assert!(report.conservation_residual.energy.abs() < 1e-11);
    }
    flow.states()
        .iter()
        .map(|state| (flow.gas().primitive(*state, || false).unwrap().pressure - 1.0).powi(2))
        .sum::<f64>()
        .sqrt()
}

#[test]
fn characteristic_outlet_lets_an_acoustic_pulse_exit_the_actual_rk_solver() {
    let characteristic = outgoing_pulse_remaining_pressure(Boundary::CharacteristicOutlet {
        pressure_pa: 1.0,
        velocity_m_per_s: 0.3,
    });
    let fixed_pressure =
        outgoing_pulse_remaining_pressure(Boundary::PressureOutlet { pressure_pa: 1.0 });
    assert!(
        fixed_pressure > 1e-4,
        "fixed-pressure reflection={fixed_pressure}"
    );
    assert!(
        characteristic < 0.02 * fixed_pressure,
        "characteristic remainder={characteristic}, fixed-pressure reflection={fixed_pressure}"
    );
}

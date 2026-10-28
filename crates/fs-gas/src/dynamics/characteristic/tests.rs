use super::*;
use crate::dynamics::tests::gas;
use crate::dynamics::{Boundaries, Boundary, Euler1d, StepControl};
use fs_qty::{Length, Time};

fn gamma_wave(pressure: f64) -> f64 {
    if pressure > 1.0 {
        (pressure - 1.0) * ((2.0 / 2.4) / (pressure + 0.4 / 2.4)).sqrt()
    } else {
        2.0 * 1.4_f64.sqrt() / 0.4 * (pressure.powf(0.4 / 2.8) - 1.0)
    }
}

#[test]
fn pressure_wave_matches_independent_gamma_hugoniot_and_rarefaction() {
    let gas = gas(0.0, 0.0);
    let inside = gas
        .primitive(gas.conserved(1.0, 0.3, 1.0).unwrap(), || false)
        .unwrap();
    for pressure in [0.05, 0.4, 1.0 - 1e-10, 1.0, 1.0 + 1e-10, 2.0, 10.0] {
        let wave = pressure_wave(&gas, inside, pressure, &mut || false).unwrap();
        let density = if pressure > 1.0 {
            (pressure + 0.4 / 2.4) / (1.0 + (0.4 / 2.4) * pressure)
        } else {
            pressure.powf(1.0 / 1.4)
        };
        let jump = gamma_wave(pressure);
        assert!(
            (wave.jump - jump).abs() < 4e-10,
            "p={pressure}, numerical={}, exact={jump}",
            wave.jump
        );
        assert!((wave.primitive.density - density).abs() < 2e-9);
        assert!((wave.primitive.velocity - (0.3 - jump)).abs() < 4e-10);
        assert!((wave.primitive.temperature - pressure / density).abs() < 3e-9);
        let exact_derivative = if pressure > 1.0 {
            ((2.0 / 2.4) / (pressure + 0.4 / 2.4)).sqrt()
                * (1.0 - 0.5 * (pressure - 1.0) / (pressure + 0.4 / 2.4))
        } else {
            pressure.powf(-(1.4 + 1.0) / 2.8) / 1.4_f64.sqrt()
        };
        assert!((wave.derivative - exact_derivative).abs() < 2e-8);
    }
}

#[test]
fn reservoir_inflow_has_original_stagnation_enthalpy_and_correct_sonic_limit() {
    let gas = gas(0.0, 0.0);
    let inside = gas.conserved(1.0, 0.0, 1.0).unwrap();
    for p0 in [1.0001, 2.0, 100.0] {
        let reservoir = gas.conserved(p0, 0.0, p0).unwrap(); // T0=1, R=1
        let critical_speed = (2.0_f64 * 1.4 / 2.4).sqrt();
        let pressure_at = |v: f64| p0 * (1.0 - v * v / 7.0).powf(3.5);
        let residual = |v: f64| v - gamma_wave(pressure_at(v));
        let mut lo = 0.0_f64;
        let mut hi = critical_speed;
        let velocity = if residual(hi) < 0.0 {
            hi
        } else {
            for _ in 0..100 {
                let middle = lo.midpoint(hi);
                if residual(middle) > 0.0 {
                    hi = middle;
                } else {
                    lo = middle;
                }
            }
            lo.midpoint(hi)
        };
        let temperature = 1.0 - velocity * velocity / 7.0;
        let pressure = pressure_at(velocity);
        let density = pressure / temperature;
        let expected_mass = -density * velocity;
        let flux = gas
            .reservoir_flux(inside, reservoir, 1.0, || false)
            .unwrap();
        assert!(
            (flux.mass - expected_mass).abs() < 2e-8 * p0,
            "p0={p0}, flux={flux:?}, expected mass={expected_mass}"
        );
        assert!((flux.momentum - (density * velocity * velocity + pressure)).abs() < 3e-8 * p0);
        assert!((flux.energy / flux.mass + gas.energy_reference() - 3.5).abs() < 2e-11);
        let mirrored = gas
            .reservoir_flux(inside, reservoir, -1.0, || false)
            .unwrap();
        assert_eq!(mirrored.mass, -flux.mass);
        assert_eq!(mirrored.energy, -flux.energy);
        assert_eq!(mirrored.momentum, flux.momentum);
    }
}

#[test]
fn outflow_samples_the_sonic_fan_and_ignores_departing_supersonic_waves() {
    let gas = gas(0.0, 0.0);
    let inside = gas.conserved(1.0, 0.0, 1.0).unwrap();
    let sink = gas.conserved(0.1, 0.0, 0.1).unwrap();
    let sonic = gas.reservoir_flux(inside, sink, 1.0, || false).unwrap();
    let speed = 2.0 * 1.4_f64.sqrt() / 2.4;
    let temperature = 4.0 / (2.4 * 2.4);
    let pressure = (2.0_f64 / 2.4).powf(2.8 / 0.4);
    let density = pressure / temperature;
    assert!((sonic.mass - density * speed).abs() < 2e-9);
    assert!((sonic.momentum - (density * speed * speed + pressure)).abs() < 2e-9);
    // The unused far-field star is colder than Tmin, but the sampled sonic
    // outlet remains inside the source interval and must still be admitted.
    let very_low = gas.conserved(1e-12, 0.0, 1e-12).unwrap();
    let same_sonic = gas.reservoir_flux(inside, very_low, 1.0, || false).unwrap();
    assert!((same_sonic.mass - sonic.mass).abs() < 2e-9);
    assert!((same_sonic.momentum - sonic.momentum).abs() < 2e-9);
    for pressure in [0.1, 2.0] {
        let q = gas.conserved(1.0, 3.0 * 1.4_f64.sqrt(), 1.0).unwrap();
        let receiver = gas.conserved(pressure, 0.0, pressure).unwrap();
        let flux = gas.reservoir_flux(q, receiver, 1.0, || false).unwrap();
        let velocity = q.momentum / q.density;
        assert_eq!(flux.mass, q.momentum);
        assert!((flux.momentum - (q.momentum * velocity + 1.0)).abs() < 2e-12);
        assert!((flux.energy - (q.energy + 1.0) * velocity).abs() < 2e-11);
    }
}

#[test]
fn supported_accelerating_inlet_does_not_require_an_unsupported_stopped_shock() {
    let base = gas(0.0, 0.0);
    let gas = FrozenGas::new(
        base.source().clone(),
        crate::dynamics::CaloricDomain::new(Temperature::new(0.001), Temperature::new(2.0)),
        || false,
    )
    .unwrap();
    let interior = gas.conserved(1.0, 0.0, 1.0).unwrap();
    let reservoir = gas.conserved(4.0, 0.0, 8.0).unwrap();
    // At p*=8 the auxiliary stopped-wave temperature is 2.2857 K; the
    // accelerating boundary has lower pressure and both actual states fit.
    assert!(
        pressure_wave(
            &gas,
            gas.primitive(interior, || false).unwrap(),
            8.0,
            &mut || false
        )
        .is_err()
    );
    let mut lo = 0.0_f64;
    let mut hi = (2.0_f64 * 1.4 * 2.0 / 2.4).sqrt();
    for _ in 0..100 {
        let v = lo.midpoint(hi);
        let p = 8.0 * (1.0 - v * v / 14.0).powf(3.5);
        if v - gamma_wave(p) > 0.0 {
            hi = v;
        } else {
            lo = v;
        }
    }
    let v = lo.midpoint(hi);
    let t = 2.0 - v * v / 7.0;
    let p = 8.0 * (t / 2.0).powf(3.5);
    let flux = gas
        .reservoir_flux(interior, reservoir, 1.0, || false)
        .unwrap();
    assert!((flux.mass + (p / t) * v).abs() < 3e-8);
    assert!((flux.energy / flux.mass + gas.energy_reference() - 7.0).abs() < 3e-11);
}

#[test]
fn reservoir_flux_retains_variable_cp_reference_and_rejects_unsupported_sources() {
    let ordinary = gas(0.4, 0.0);
    let shifted = gas(0.4, 1e18);
    for (rho, velocity, pressure, reservoir_pressure) in [
        (1.0, 0.3, 1.0, 0.7),
        (1.0, 0.0, 1.0, 2.0),
        (1.0, -0.2, 1.0, 1.0),
    ] {
        let q = ordinary.conserved(rho, velocity, pressure).unwrap();
        let reservoir = ordinary
            .conserved(reservoir_pressure / 2.0, 0.0, reservoir_pressure)
            .unwrap();
        let a = ordinary
            .reservoir_flux(q, reservoir, 1.0, || false)
            .unwrap();
        let b = shifted.reservoir_flux(q, reservoir, 1.0, || false).unwrap();
        assert_eq!(a, b);
        if a.mass < 0.0 {
            // Manufactured source h(T)=3.5T+0.2T², R=1; formation=0.
            assert!((a.energy / a.mass + ordinary.energy_reference() - 7.8).abs() < 3e-11);
        }
    }
    let q = ordinary.conserved(1.0, 0.0, 1.0).unwrap();
    assert_eq!(
        ordinary.reservoir_flux(q, q, 1.0, || true),
        Err(Error::Cancelled)
    );
    let moving = ordinary.conserved(1.0, 0.1, 1.0).unwrap();
    assert!(matches!(
        ordinary.reservoir_flux(q, moving, 1.0, || false),
        Err(Error::Invalid(_))
    ));
    assert!(matches!(
        ordinary.reservoir_flux(q, q, 0.0, || false),
        Err(Error::Invalid(_))
    ));
    let cold = ordinary.conserved(2000.0, 0.0, 2.0).unwrap();
    assert!(matches!(
        ordinary.reservoir_flux(q, cold, 1.0, || false),
        Err(Error::CaloricPath(_))
    ));
}

#[test]
fn driven_euler_boundaries_publish_actual_end_exchange_and_update_atomically() {
    let gas = gas(0.3, 0.0);
    let ambient = gas.conserved(1.0, 0.0, 1.0).unwrap();
    let feed = gas.conserved(1.0, 0.0, 1.2).unwrap();
    let boundaries = Boundaries::Ends {
        left: Boundary::Reservoir(feed),
        right: Boundary::Reservoir(ambient),
    };
    let mut flow = Euler1d::new(
        gas.clone(),
        Length::new(0.025),
        0.1,
        vec![ambient; 40],
        boundaries,
        || false,
    )
    .unwrap();
    let initial = flow.totals();
    let mut exchange = [0.0; 2];
    for _ in 0..60 {
        let report = flow
            .step(Time::new(0.001), StepControl::default(), || false)
            .unwrap();
        exchange[0] += report.left_boundary_change.mass + report.right_boundary_change.mass;
        exchange[1] += report.left_boundary_change.energy + report.right_boundary_change.energy;
        assert!(report.conservation_residual.mass.abs() < 2e-14);
        assert!(report.conservation_residual.energy.abs() < 2e-13);
    }
    assert!(exchange[0] > 0.0 && exchange[1] > 0.0);
    assert!((flow.totals().mass - initial.mass - exchange[0]).abs() < 2e-13);
    assert!((flow.totals().energy - initial.energy - exchange[1]).abs() < 2e-12);
    let before = flow.states().to_vec();
    let time = flow.time();
    let bad = gas.conserved(1.0, 0.1, 1.0).unwrap();
    assert!(
        flow.set_boundaries(
            Boundaries::Ends {
                left: Boundary::Reservoir(bad),
                right: Boundary::Reflecting
            },
            || false
        )
        .is_err()
    );
    assert_eq!(flow.boundaries(), boundaries);
    assert_eq!(
        flow.set_boundaries(Boundaries::Periodic, || true),
        Err(Error::Cancelled)
    );
    assert_eq!(flow.boundaries(), boundaries);
    assert_eq!(flow.states(), before);
    assert_eq!(flow.time(), time);
    flow.set_boundaries(
        Boundaries::Ends {
            left: Boundary::Reflecting,
            right: Boundary::Reflecting,
        },
        || false,
    )
    .unwrap();
    let report = flow
        .step(Time::new(0.001), StepControl::default(), || false)
        .unwrap();
    assert_eq!(report.boundary_change.mass, 0.0);
    assert_eq!(report.boundary_change.energy, 0.0);
}

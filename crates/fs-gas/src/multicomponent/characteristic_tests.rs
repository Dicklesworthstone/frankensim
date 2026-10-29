//! Independent donor-composition and thermodynamic checks for compact ports.
use super::tests::{gas, source};
use super::{Conserved, Error, MixtureGas};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

pub(super) fn caloric_gas(variable: bool, formation: f64) -> MixtureGas {
    let initial = gas(vec![
        source("A", R, 3.5, if variable { 0.25 } else { 0.0 }, formation),
        source(
            "B",
            2.0 * R,
            5.0,
            if variable { 0.7 } else { 0.0 },
            -2.0 * formation,
        ),
    ]);
    let mut domain = initial.domain();
    // Resolve the weak-wave oracle independently of the general test fixture's
    // looser 1e-10 J/kg absolute caloric inversion budget.
    domain.absolute_energy_tolerance = 1e-14;
    domain.relative_energy_tolerance = 1e-14;
    MixtureGas::new(initial.source().clone(), domain, || false).unwrap()
}

fn flux(gas: &MixtureGas, inside: &Conserved, reservoir: &Conserved, normal: f64) -> [f64; 4] {
    let mut values = [f64::NAN; 4];
    let receipt = gas
        .reservoir_flux(
            inside.as_cell(),
            reservoir.as_cell(),
            normal,
            &mut values,
            || false,
        )
        .unwrap();
    assert!(receipt.signal_speed.is_finite() && receipt.signal_speed > 0.0);
    assert!(!receipt.used_hll_fallback);
    assert!(values.iter().all(|v| v.is_finite()));
    values
}
fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "actual={actual:.15e}, expected={expected:.15e}, scale={scale:.9e}, tolerance={tolerance:.3e}"
    );
}
fn gamma_wave(pressure: f64) -> f64 {
    // The interior trace is pure A: rho=p=1, gamma=1.4.
    if pressure > 1.0 {
        (pressure - 1.0) * ((2.0 / 2.4) / (pressure + 0.4 / 2.4)).sqrt()
    } else {
        2.0 * 1.4_f64.sqrt() / 0.4 * (pressure.powf(0.4 / 2.8) - 1.0)
    }
}

#[test]
fn different_compositions_match_independent_inlet_and_sonic_mass_enthalpy() {
    let gas = caloric_gas(false, 0.0);
    let inside = gas.conserved(&[1.0, 0.0], 0.0, 1.0).unwrap();
    // Donor B has R=1/2, Cp=5/2, gamma=5/4 and T0=2. The duct wave
    // instead uses species A's gamma=7/5; freezing either composition into
    // the other gives a different characteristic match and mass flux.
    for p0 in [1.0001, 2.0, 100.0] {
        let reservoir = gas.conserved(&[0.0, p0], 0.0, p0).unwrap();
        let critical = (2.0_f64 * 1.25 * 0.5 * 2.0 / 2.25).sqrt();
        let pressure_at = |v: f64| p0 * (1.0 - v * v / 10.0).powf(5.0);
        let residual = |v: f64| v - gamma_wave(pressure_at(v));
        let (mut lo, mut hi) = (0.0_f64, critical);
        let speed = if residual(hi) < 0.0 {
            hi
        } else {
            for _ in 0..100 {
                let v = lo.midpoint(hi);
                if residual(v) > 0.0 {
                    hi = v;
                } else {
                    lo = v;
                }
            }
            lo.midpoint(hi)
        };
        let temperature = 2.0 - speed * speed / 5.0;
        let pressure = pressure_at(speed);
        let density = pressure / (0.5 * temperature);
        let expected_mass = -density * speed;
        let actual = flux(&gas, &inside, &reservoir, 1.0);
        assert_eq!(
            actual[0], 0.0,
            "non-donor species A entered from a pure-B reservoir"
        );
        close(actual[1], expected_mass, p0, 2e-10);
        close(actual[2], density * speed * speed + pressure, p0, 2e-10);
        // Stored energy removes e_B(Tmin)=Cv_B*Tmin. Total stagnation
        // enthalpy is otherwise exactly Cp_B*T0, not interior A's enthalpy.
        close(actual[3] / actual[1], 2.5 * 2.0 - 2.0 * 0.001, 5.0, 2e-12);
        if p0 == 100.0 {
            close(speed, critical, critical, 1e-14);
        }
    }
}

#[test]
fn sonic_outflow_carries_interior_species_and_ignores_unused_cold_tail() {
    let gas = caloric_gas(false, 0.0);
    for densities in [[1.0_f64, 0.0], [0.3, 0.7]] {
        let r = densities[0] + 0.5 * densities[1];
        let cp = 3.5 * densities[0] + 2.5 * densities[1];
        let cv = cp - r;
        let gamma = cp / cv;
        let factor = 2.0 / (gamma + 1.0);
        let speed = factor * gamma.sqrt();
        let temperature = factor * factor / r;
        let pressure = factor.powf(2.0 * gamma / (gamma - 1.0));
        let density = pressure / (r * temperature);
        let expected_mass = density * speed;
        let inside = gas.conserved(&densities, 0.0, 1.0).unwrap();
        for backpressure in [0.02, 1e-12] {
            let sink = gas
                .conserved(&[0.0, backpressure], 0.0, backpressure)
                .unwrap();
            let actual = flux(&gas, &inside, &sink, 1.0);
            for k in 0..2 {
                close(actual[k], densities[k] * expected_mass, 1.0, 3e-10);
            }
            close(actual[2], density * speed * speed + pressure, 1.0, 3e-10);
            let enthalpy = cp * temperature - cv * 0.001 + 0.5 * speed * speed;
            close(actual[3], expected_mass * enthalpy, 4.0, 3e-10);
            if densities[1] == 0.0 {
                assert_eq!(actual[1], 0.0);
            }
        }
    }
}

#[test]
fn weak_mixture_inflow_has_the_branch_acoustic_impedance() {
    let gas = caloric_gas(false, 0.0);
    let inside = gas.conserved(&[1.0, 0.0], 0.0, 1.0).unwrap();
    let sound = 1.4_f64.sqrt();
    for excess in [1e-5, 1e-6, 1e-7] {
        let reservoir = gas
            .conserved(&[0.0, 1.0 + excess], 0.0, 1.0 + excess)
            .unwrap();
        let actual = flux(&gas, &inside, &reservoir, 1.0);
        // p'=rho_i*c_i*v at leading order. The stationary reservoir adds
        // no second tube impedance; donor rho0=1 at the expansion point.
        close(-actual[1] * sound / excess, 1.0, 1.0, 3e-5);
        assert_eq!(actual[0], 0.0);
    }
    let equilibrium = gas.conserved(&[0.0, 1.0], 0.0, 1.0).unwrap();
    let stationary = flux(&gas, &inside, &equilibrium, 1.0);
    for k in [0, 1, 3] {
        close(stationary[k], 0.0, 1.0, 2e-12);
    }
    close(stationary[2], 1.0, 1.0, 2e-12);
}

#[test]
fn variable_cp_ports_preserve_donor_enthalpy_reference_and_orientation() {
    let ordinary = caloric_gas(true, 0.0);
    for (inside_d, velocity, pressure, reservoir_d, reservoir_pressure) in [
        ([0.8, 0.2], 0.0, 1.0, [0.3, 0.7], 2.0),
        ([0.2, 0.8], 0.3, 1.0, [0.9, 0.1], 0.6),
    ] {
        let inside = ordinary.conserved(&inside_d, velocity, pressure).unwrap();
        let reservoir = ordinary
            .conserved(&reservoir_d, 0.0, reservoir_pressure)
            .unwrap();
        let expected = flux(&ordinary, &inside, &reservoir, 1.0);
        let mass = expected[0] + expected[1];
        let donor = if mass < 0.0 { reservoir_d } else { inside_d };
        for k in 0..2 {
            close(expected[k], mass * donor[k], mass.abs().max(1.0), 2e-13);
        }
        if mass < 0.0 {
            let r = donor[0] + 0.5 * donor[1];
            let a = 3.5 * donor[0] + 2.5 * donor[1];
            let b = 0.25 * donor[0] + 0.35 * donor[1];
            let t0 = reservoir_pressure / r;
            let h0 = a * t0 + 0.5 * b * t0 * t0;
            let reference = (a - r) * 0.001 + 0.5 * b * 0.001_f64.powi(2);
            close(expected[3] / mass, h0 - reference, h0, 2e-12);
        }
        for reference in [-1e18, 0.0, 1e18] {
            let shifted = caloric_gas(true, reference);
            for normal in [-1.0, 1.0] {
                let inside = shifted
                    .conserved(&inside_d, normal * velocity, pressure)
                    .unwrap();
                let reservoir = shifted
                    .conserved(&reservoir_d, 0.0, reservoir_pressure)
                    .unwrap();
                let actual = flux(&shifted, &inside, &reservoir, normal);
                assert_eq!(
                    actual,
                    [
                        normal * expected[0],
                        normal * expected[1],
                        expected[2],
                        normal * expected[3]
                    ]
                );
            }
        }
    }
}

#[test]
fn port_refusal_and_late_cancellation_leave_the_entire_species_buffer_unchanged() {
    let gas = caloric_gas(true, 0.0);
    let inside = gas.conserved(&[0.8, 0.2], 0.0, 1.0).unwrap();
    let reservoir = gas.conserved(&[0.3, 0.7], 0.0, 2.0).unwrap();
    let mut probe = [0.0; 4];
    let mut count = 0;
    gas.reservoir_flux(
        inside.as_cell(),
        reservoir.as_cell(),
        1.0,
        &mut probe,
        || {
            count += 1;
            false
        },
    )
    .unwrap();
    assert!(count > 10);
    for stop in [1, count / 2, count] {
        let mut output = [7.0; 4];
        let mut visits = 0;
        let result = gas.reservoir_flux(
            inside.as_cell(),
            reservoir.as_cell(),
            1.0,
            &mut output,
            || {
                visits += 1;
                visits >= stop
            },
        );
        assert!(matches!(result, Err(Error::Cancelled)), "stop={stop}");
        assert_eq!(output, [7.0; 4]);
    }
    let moving = gas.conserved(&[0.3, 0.7], 0.1, 2.0).unwrap();
    for (outside, normal) in [(&moving, 1.0), (&reservoir, 0.0)] {
        let mut output = [7.0; 4];
        assert!(
            gas.reservoir_flux(
                inside.as_cell(),
                outside.as_cell(),
                normal,
                &mut output,
                || false
            )
            .is_err()
        );
        assert_eq!(output, [7.0; 4]);
    }
    let mut short = [7.0; 3];
    assert!(
        gas.reservoir_flux(
            inside.as_cell(),
            reservoir.as_cell(),
            1.0,
            &mut short,
            || false
        )
        .is_err()
    );
    assert_eq!(short, [7.0; 3]);
}

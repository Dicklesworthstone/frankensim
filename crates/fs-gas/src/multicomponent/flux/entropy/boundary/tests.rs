//! Independent moving-exterior entropy effort checks in the source energy chart.

use super::*;
use crate::multicomponent::tests::{gas, source};

// Cp/R = a+b*T; molar mass = mass*R, reference pressure = 1.
const ROWS: [[f64; 3]; 3] = [[3.5, 0.2, 1.0], [4.0, 0.4, 2.0], [3.5, 0.0, 1.0]];
const T_MIN: f64 = 0.001;
const T_REF: f64 = 1.3;

fn density_basis() -> [f64; MAX_SPECIES] {
    let mut basis = [1.0; MAX_SPECIES];
    basis[..3].copy_from_slice(&[0.8, 0.2, 1.0]);
    basis
}

fn alpha_per_mass(k: usize) -> f64 {
    let [a, b, mass] = ROWS[k];
    ((a - 1.0) * T_REF.ln() + b * T_REF - (density_basis()[k] / mass).ln() - 1.0) / mass
}

fn fixture(formation: f64, velocity: f64) -> (MixtureGas, crate::multicomponent::Conserved) {
    let gas = gas(vec![
        source("A", R, 3.5, 0.2, formation),
        source("B", 2.0 * R, 4.0, 0.4, -2.0 * formation),
        source("C", R, 3.5, 0.0, 3.0 * formation),
    ]);
    let rho = [0.6, 1.1, 0.0];
    let q = gas.conserved(&rho, velocity, 3.2 * (rho[0] + 0.5 * rho[1])).unwrap();
    (gas, q)
}

fn near(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 2e-10 * expected.abs().max(1.0),
        "actual={actual:e}, expected={expected:e}");
}

#[test]
fn g0_moving_exterior_physical_flux_recovers_advected_source_entropy_and_source_gauge() {
    for velocity in [2.3, -1.7] {
        let mut original = None;
        for formation in [0.0, 1e18] {
            let (gas, q) = fixture(formation, velocity);
            let p = gas.primitive(q.as_cell(), || false).unwrap();
            let effort = ExteriorEntropyEffort::new(
                &gas, q.as_cell(), p, T_REF, &density_basis(), &mut || false,
            ).unwrap();
            let mut flux = [0.0; 5];
            gas.entropy_stable_flux(q.as_cell(), q.as_cell(), &mut flux, || false).unwrap();
            let source_entropy: f64 = (0..2).map(|k| {
                let [a, b, mass] = ROWS[k];
                let rho = q.partial_densities[k];
                rho / mass * ((a - 1.0) * p.temperature.ln() + b * p.temperature
                    - (rho / mass).ln())
            }).sum();
            for signed_measure in [0.17, -0.17] {
                let transfer = flux.map(|value| signed_measure * value);
                let input = effort.entropy_input(&transfer, signed_measure, &mut || false)
                    .unwrap().to_f64();
                let restored: f64 = (0..2).map(|k| alpha_per_mass(k) * transfer[k]).sum();
                near(input + restored, signed_measure * p.velocity * source_entropy);
            }
            let input = effort.entropy_input(&flux, 1.0, &mut || false).unwrap().to_f64();
            if let Some(expected) = original {
                assert_eq!(input, expected);
            } else {
                original = Some(input);
            }
        }
    }
}

#[test]
fn g0_arbitrary_transfer_retains_momentum_kinetic_and_entropy_potential_terms() {
    let (gas, q) = fixture(0.0, 2.3);
    let p = gas.primitive(q.as_cell(), || false).unwrap();
    let effort = ExteriorEntropyEffort::new(
        &gas, q.as_cell(), p, T_REF, &density_basis(), &mut || false,
    ).unwrap();
    let transfer = [0.012, -0.004, 0.0, 0.03, 0.07];
    let measure = 0.05;
    let mut expected = transfer[4] / p.temperature - p.velocity * transfer[3] / p.temperature;
    let mut stationary_terms = transfer[4] / p.temperature;
    for k in 0..2 {
        let [a, b, mass] = ROWS[k];
        let t = p.temperature;
        let energy = ((a - 1.0) * (t - T_MIN)
            + 0.5 * b * (t * t - T_MIN * T_MIN)) / mass;
        let entropy = ((a - 1.0) * t.ln() + b * t
            - (q.partial_densities[k] / mass).ln()) / mass;
        let chemical = energy / t + 1.0 / mass - entropy + alpha_per_mass(k);
        expected -= (chemical - 0.5 * p.velocity * p.velocity / t) * transfer[k];
        stationary_terms -= chemical * transfer[k];
        expected += measure * p.velocity * q.partial_densities[k] / mass;
    }
    let input = effort.entropy_input(&transfer, measure, &mut || false).unwrap().to_f64();
    near(input, expected);
    assert!((input - stationary_terms).abs() > 0.05);
    let reversed = transfer.map(|value| -value);
    near(effort.entropy_input(&reversed, -measure, &mut || false).unwrap().to_f64(), -input);
}

#[test]
fn g4_absent_activity_is_checked_before_stage_cancellation_and_bad_axes_refuse() {
    let (gas, q) = fixture(0.0, 0.7);
    let p = gas.primitive(q.as_cell(), || false).unwrap();
    let effort = ExteriorEntropyEffort::new(
        &gas, q.as_cell(), p, T_REF, &density_basis(), &mut || false,
    ).unwrap();
    for transfer in [1e-100, -1e-100] {
        let stage = [0.0, 0.0, transfer, 0.0, 0.0];
        assert!(matches!(effort.admit_flux(&stage, &mut || false), Err(Error::EntropyDomain(_))));
        assert!(matches!(effort.entropy_input(&stage, 0.0, &mut || false), Err(Error::EntropyDomain(_))));
    }
    effort.admit_flux(&[0.0; 5], &mut || false).unwrap();
    assert!(matches!(effort.admit_flux(&[0.0; 4], &mut || false), Err(Error::EntropyDomain(_))));
    assert!(matches!(
        effort.entropy_input(&[0.0; 5], f64::INFINITY, &mut || false),
        Err(Error::EntropyDomain(_)),
    ));
    assert!(matches!(
        effort.admit_flux(&[0.0, 0.0, 0.0, f64::NAN, 0.0], &mut || false),
        Err(Error::EntropyDomain(_)),
    ));
    assert_eq!(effort.admit_flux(&[0.0; 5], &mut || true), Err(Error::Cancelled));
    assert_eq!(
        effort.entropy_input(&[0.0; 5], 0.0, &mut || true).unwrap_err(),
        Error::Cancelled,
    );
}

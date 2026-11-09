use super::*;
use crate::multicomponent::tests::{gas, source};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn model(formation: f64) -> MixtureGas {
    gas(vec![
        source("A", R, 3.5, 0.1, formation),
        source("B", 2.0 * R, 4.5, 0.2, -2.0 * formation),
        source("absent-C", R, 3.5, 0.0, 7.0 * formation),
    ])
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_mixture_storage_mixes_species_enthalpy_and_full_incoming_kinetic_energy() {
    let mut reference = None;
    for formation in [0.0, -1e18, 1e18] {
        let gas = model(formation);
        let volume = 0.4;
        let tank = Plenum::new(volume, gas.conserved(&[1.0, 0.0, 0.0], 0.0, 2.0).unwrap()).unwrap();
        let donor = gas.conserved(&[0.0, 1.0, 0.0], 3.0, 1.5).unwrap();
        let masses = [0.0, 0.1, 0.0];
        let energy = masses[1] * (donor.energy + 1.5);
        let candidate = tank
            .exchange_candidate(&gas, &masses, energy, || false)
            .unwrap();
        let mut output = tank.clone();
        let pointer = output.state().partial_densities.as_ptr();
        tank.exchange_candidate_into(&gas, &masses, energy, &mut output, || false)
            .unwrap();
        assert_eq!(output, candidate);
        assert_eq!(output.state().partial_densities.as_ptr(), pointer);
        assert_eq!(candidate.state().partial_densities[2], 0.0);
        assert_eq!(candidate.state().momentum, 0.0);
        close(candidate.state().partial_densities[0] * volume, 0.4, 2e-15);
        close(candidate.state().partial_densities[1] * volume, 0.1, 2e-15);
        let total_energy = volume * tank.state().energy + energy;
        close(volume * candidate.state().energy, total_energy, 2e-15);
        // Independent integral of sum(m_i Cv_i(T)), with R_A=1, R_B=1/2.
        let a = 0.4 * 2.5 + 0.1 * 0.5 * 3.5;
        let b = 0.5 * (0.4 * 0.1 + 0.1 * 0.5 * 0.2);
        let c = total_energy + a * 0.001 + b * 0.001_f64.powi(2);
        let temperature = 2.0 * c / (a + (a * a + 4.0 * b * c).sqrt());
        let p = gas.primitive(candidate.state(), || false).unwrap();
        close(p.temperature, temperature, 2e-10);
        close(p.pressure, (0.4 + 0.1 * 0.5) * temperature / volume, 2e-10);
        if let Some(previous) = &reference {
            assert_eq!(&candidate, previous);
        } else {
            reference = Some(candidate.clone());
        }
        assert_eq!(
            tank.exchange_candidate(&gas, &[0.0; 3], 0.0, || false)
                .unwrap(),
            tank
        );
    }
}

#[test]
fn g4_mixture_output_survives_species_exhaustion_shape_source_and_late_cancel() {
    let gas = model(0.0);
    let tank = Plenum::new(0.4, gas.conserved(&[0.25, 0.75, 0.0], 0.0, 1.0).unwrap()).unwrap();
    let initial = gas.primitive(tank.state(), || false).unwrap();
    let pure = gas
        .conserved(&[1.0, 0.0, 0.0], 0.0, initial.temperature)
        .unwrap();
    let exhausted = tank
        .exchange_candidate(&gas, &[-0.1, 0.0, 0.0], -0.1 * pure.energy, || false)
        .unwrap();
    assert_eq!(exhausted.state().partial_densities[0], 0.0);
    assert_eq!(exhausted.state().partial_densities[2], 0.0);
    close(
        gas.primitive(exhausted.state(), || false)
            .unwrap()
            .temperature,
        initial.temperature,
        2e-10,
    );

    let mut output = Plenum::new(0.7, gas.conserved(&[0.0, 1.0, 0.0], 0.0, 2.0).unwrap()).unwrap();
    let before = output.clone();
    for (masses, energy) in [
        (vec![0.0; 2], 0.0),
        (vec![0.0; 4], 0.0),
        (vec![-0.11, 0.0, 0.0], 0.0),
        (vec![f64::NAN, 0.0, 0.0], 0.0),
        (vec![0.0; 3], -100.0),
        (vec![0.0; 3], 1e8),
        (vec![0.0; 3], f64::INFINITY),
    ] {
        assert!(
            tank.exchange_candidate_into(&gas, &masses, energy, &mut output, || false)
                .is_err()
        );
        assert_eq!(output, before);
    }
    let mut short = Plenum::new(
        0.4,
        Conserved {
            partial_densities: vec![1.0],
            momentum: 0.0,
            energy: 1.0,
        },
    )
    .unwrap();
    let short_before = short.clone();
    assert!(
        tank.exchange_candidate_into(&gas, &[0.0; 3], 0.0, &mut short, || false)
            .is_err()
    );
    assert_eq!(short, short_before);
    let mut polls = 0;
    tank.exchange_candidate_into(&gas, &[0.01, 0.02, 0.0], 0.2, &mut output, || {
        polls += 1;
        false
    })
    .unwrap();
    let accepted = output.clone();
    output = before.clone();
    let pointer = output.state().partial_densities.as_ptr();
    let mut calls = 0;
    assert!(matches!(
        tank.exchange_candidate_into(&gas, &[0.01, 0.02, 0.0], 0.2, &mut output, || {
            calls += 1;
            calls == polls
        }),
        Err(Error::Cancelled)
    ));
    assert_eq!(output, before);
    assert_eq!(output.state().partial_densities.as_ptr(), pointer);
    tank.exchange_candidate_into(&gas, &[0.01, 0.02, 0.0], 0.2, &mut output, || false)
        .unwrap();
    assert_eq!(output, accepted);
}

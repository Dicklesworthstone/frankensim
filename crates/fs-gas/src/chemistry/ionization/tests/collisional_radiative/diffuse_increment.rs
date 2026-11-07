use super::*;
use crate::chemistry::collisional_radiative::solve::{DiffuseIncrement, Problem};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

const CV: f64 = 1.5 * R;

#[test]
fn g0_opaque_recycling_retains_independent_net_energy() {
    let (initial, source) = fixture(0.5, 10000.0, 0.0, 0.0, 0.0);
    let controls = controls(1e-7);
    let problem = Problem::new(&initial, &source, &[], controls, &mut || false).unwrap();
    let point = problem.initial_point();
    let chi = source.source().channel().ionization_energy_j_mol();
    // Many self-recycling events can coexist with a finite, far smaller loss
    // through the boundary. Manufacture the independently transported net
    // moments, below the floating-point resolution of either gross account.
    let extent = -1e-3;
    let kinetic_per_mol = 0.6 * CV * 10000.0;
    let electron_energy = extent * kinetic_per_mol;
    let net_energy = chi * extent + electron_energy;
    let capture = 1e15;
    let emitted = chi * capture + kinetic_per_mol * capture;
    let increment = DiffuseIncrement {
        net_extent_mol: extent,
        net_electron_energy_j: electron_energy,
        net_gas_energy_j: net_energy,
        absorbed_mol: capture + extent,
        captured_mol: capture,
        absorbed_energy_j: emitted + net_energy,
        captured_energy_j: kinetic_per_mol * capture,
    };
    assert_eq!(increment.absorbed_mol, increment.captured_mol);
    assert_eq!(increment.absorbed_energy_j, emitted);
    let next = problem.apply_diffuse(point, increment).unwrap();
    let expected_temperature = (CV * 0.5 * 10000.0 + electron_energy) / (CV * 0.499);
    assert!((next.electron_temperature() - expected_temperature).abs() < 2e-11);
    assert_eq!(
        next.heavy_temperature(),
        initial.heavy_temperature().value()
    );
    assert_eq!(next.net_radiation_energy_j(), net_energy);
    assert_eq!(problem.active_amounts(next), [0.501, 0.499, 0.499]);
    let mut balance = CollisionalRadiativeBalance::default();
    let state = problem.publish(next, &mut balance, &mut || false).unwrap();
    assert_eq!(state.amounts()[1], state.amounts()[2]);
    assert!((balance.electron_sensible_energy_change_j - electron_energy).abs() < 1e-9);
    assert!(
        (balance.chemical_storage_change_j + balance.electron_sensible_energy_change_j
            - net_energy)
            .abs()
            < 1e-9
    );
    assert!(balance.energy_residual_j.abs() <= controls.energy_tolerance_j);
    assert_eq!(balance.net_radiation_energy_j, net_energy);
    assert!(balance.heavy_sensible_energy_residual_j.abs() <= controls.energy_tolerance_j);
}

#[test]
fn g0_diffuse_first_electrons_and_trial_refusal_preserve_state() {
    let (initial, source) = fixture(0.0, 1000.0, 0.0, 0.0, 0.0);
    let problem = Problem::new(&initial, &source, &[], controls(1e-7), &mut || false).unwrap();
    let point = problem.initial_point();
    let extent = 0.01;
    let electron_energy = CV * extent * 16000.0;
    let absorbed = source.source().channel().ionization_energy_j_mol() * extent + electron_energy;
    let increment = DiffuseIncrement {
        net_extent_mol: extent,
        net_electron_energy_j: electron_energy,
        net_gas_energy_j: absorbed,
        absorbed_mol: extent,
        absorbed_energy_j: absorbed,
        ..DiffuseIncrement::default()
    };
    let next = problem.apply_diffuse(point, increment).unwrap();
    assert!((next.electron_temperature() - 16000.0).abs() < 1e-11);
    assert_eq!(problem.active_amounts(next), [0.99, 0.01, 0.01]);
    assert!(matches!(
        problem.apply_diffuse(
            point,
            DiffuseIncrement {
                net_extent_mol: -1e-6,
                ..DiffuseIncrement::default()
            }
        ),
        Err(CollisionalRadiativeError::OutsideDomain)
    ));
    assert!(matches!(
        problem.apply_diffuse(
            point,
            DiffuseIncrement {
                net_gas_energy_j: absorbed + 1.0,
                ..increment
            }
        ),
        Err(CollisionalRadiativeError::Invalid(_))
    ));
    assert_eq!(problem.active_amounts(point), [1.0, 0.0, 0.0]);
    assert_eq!(point.net_radiation_energy_j(), 0.0);
}

#[test]
fn g0_collision_subflow_excludes_both_photon_processes() {
    let (initial, source) = fixture(0.1, 10000.0, 0.0, 0.0, 0.0);
    let photons = field(20.0 * EV, 1000.0);
    let problem = Problem::new(&initial, &source, &photons, controls(1e-7), &mut || false).unwrap();
    let point = problem.initial_point();
    let mut balance = CollisionalRadiativeBalance::default();
    let rates = problem
        .collisional_rates(point, &mut balance, &mut || false)
        .unwrap();
    assert!(rates.inactive(point));
    let next = problem
        .collisional_midpoint(point, rates, 1.0, &mut balance, &mut || false)
        .unwrap();
    let state = problem.publish(next, &mut balance, &mut || false).unwrap();
    assert_eq!(state, initial);
    assert_eq!(balance.absorbed_photons, 0.0);
    assert_eq!(balance.emitted_photons, 0.0);
    assert_eq!(balance.maximum_relative_quadrature_error, 0.0);
    let illuminated = problem.rates(point, &mut balance, &mut || false).unwrap();
    assert!(!illuminated.inactive(point));
}

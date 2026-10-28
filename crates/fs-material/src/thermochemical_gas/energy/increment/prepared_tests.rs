use super::*;
use crate::thermochemical_gas::tests::{model, point};

fn controls() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 }
}

#[test]
fn prepared_reuse_matches_original_complete_pressure_and_density_solutions() {
    let gas = model(vec![0.75, 0.25], 0.002);
    let branch = gas.prepare_caloric_branch(controls(), || false).unwrap();
    assert!(core::ptr::eq(branch.model(), &gas));
    assert_eq!(branch.controls(), controls());
    for temperature in [300.0, 500.0, 1000.0] {
        let context = point(temperature, 100000.0).with("material-lot", 42.0).unwrap();
        let state = gas.evaluate(&context).unwrap();
        let density = state.state().density;
        assert!(branch.minimum_cp_j_kg_k() <= state.state().specific_heat_cp);
        assert!(branch.minimum_cv_j_kg_k() <= state.thermodynamics().properties().mass_specific_cv().value());
        for change in [-10000.0, 0.0, 10000.0] {
            assert_eq!(branch.reconstruct_density_internal_energy_increment(
                &context, density, change, || false).unwrap(),
                gas.reconstruct_density_internal_energy_increment(
                    &context, density, change, controls(), || false).unwrap());
            assert_eq!(branch.reconstruct_pressure_enthalpy_increment(
                &context, change, || false).unwrap(),
                gas.reconstruct_pressure_enthalpy_increment(
                    &context, change, controls(), || false).unwrap());
        }
    }
}

#[test]
fn prepared_admission_is_amortized_without_weakening_cancellation_or_domains() {
    let gas = model(vec![0.75, 0.25], 0.002);
    let mut admission_polls = 0;
    let branch = gas.prepare_caloric_branch(controls(), || { admission_polls += 1; false }).unwrap();
    let mut cached_polls = 0;
    let context = point(400.0, 100000.0);
    let expected = branch.reconstruct_pressure_enthalpy_increment(&context, 10000.0,
        || { cached_polls += 1; false }).unwrap();
    let mut ordinary_polls = 0;
    assert_eq!(expected, gas.reconstruct_pressure_enthalpy_increment(&context, 10000.0, controls(),
        || { ordinary_polls += 1; false }).unwrap());
    assert_eq!(ordinary_polls, admission_polls + cached_polls);
    assert!(matches!(gas.prepare_caloric_branch(controls(), || true), Err(GasEnergyError::Cancelled)));
    for stop in 1..=cached_polls {
        let mut polls = 0;
        assert!(matches!(branch.reconstruct_pressure_enthalpy_increment(&context, 10000.0,
            || { polls += 1; polls == stop }), Err(GasEnergyError::Cancelled)));
    }
    assert_eq!(expected, branch.reconstruct_pressure_enthalpy_increment(&context, 10000.0, || false).unwrap());
    assert!(branch.reconstruct_pressure_enthalpy_increment(&point(1900.0, 100000.0), 0.0, || false).is_err());
    assert!(gas.prepare_caloric_branch(GasEnergySolve { maximum_temperature_k: 2100.0, ..controls() }, || false).is_err());
}

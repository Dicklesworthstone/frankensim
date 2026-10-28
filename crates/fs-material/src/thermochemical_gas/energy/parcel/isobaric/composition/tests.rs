use super::*;
use crate::thermochemical_gas::tests::{close, point, species_in_domain, transport_in_domain};
use fs_thermochem::{Composition, CompositionBasis, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;

fn controls() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1500.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 }
}
fn catalog(product_minimum: f64) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(vec![
        species_in_domain("A", 0.028, 3.5, 0.001, 0.0, [200.0, 2000.0]),
        species_in_domain("B", 0.056, 6.5, 0.002, -300000.0, [product_minimum, 2000.0]),
    ], Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5]).unwrap()).unwrap(), vec![
        transport_in_domain("A", 1.7e-5, 0.025, [200.0, 2000.0], [1000.0, 1e7]),
        transport_in_domain("B", 2e-5, 0.15, [product_minimum, 2000.0], [1000.0, 1e7]),
    ], None).unwrap()
}
fn composition(source: &ThermochemicalGasModel, fractions: [f64; 2]) -> ThermochemicalGasModel {
    source.with_composition(Composition::new(CompositionBasis::MoleFraction, fractions.to_vec()).unwrap()).unwrap()
}
fn h_a(t: f64) -> f64 { (R / 0.028) * (3.5 * t + 0.0005 * t * t) }
fn h_b(t: f64) -> f64 { (R / 0.056) * (6.5 * t + 0.001 * t * t) - 300000.0 }
fn original(source: &ThermochemicalGasModel) -> ClosedGasParcel {
    ClosedGasParcel::new(composition(source, [1.0, 0.0]),
        point(500.0, 100000.0).with("material-lot", 42.0).unwrap(), 0.01).unwrap()
}

#[test]
fn isobaric_composition_heat_uses_formation_enthalpy_and_new_molecular_amount() {
    let source = catalog(200.0);
    for target in [700.0, 1100.0] {
        let mut gas = original(&source);
        let before = gas.clone();
        let heat = gas.mass_kg() * (h_b(target) - h_a(500.0));
        let step = gas.replace_model_with_heat_at_constant_pressure(
            composition(&source, [0.0, 1.0]), heat, controls(), || false,
        ).unwrap();
        close(gas.state().state().temperature, target);
        assert_eq!(gas.state().state().pressure, before.state().state().pressure);
        assert_eq!(gas.mass_kg(), before.mass_kg());
        close(gas.volume_m3(), 0.01 * 0.5 * target / 500.0);
        let work = 100000.0 * (gas.volume_m3() - before.volume_m3());
        assert_eq!(step.energy.work_by_gas_j, work);
        close(step.enthalpy_change_j, heat);
        close(step.energy.stored_energy_change_j, heat - work);
        assert!(step.energy.energy_residual_j.abs() <= gas.mass_kg() * controls().energy_tolerance_j_kg);
        assert!(step.enthalpy_residual_j.abs() <= gas.mass_kg() * controls().energy_tolerance_j_kg);
        assert_eq!(gas.point().axes().get("material-lot"), Some(&42.0));
        assert_eq!(gas.state().components().len(), 1);
        gas.state().state().validate_coherence().unwrap();
    }
}

#[test]
fn adiabatic_two_to_one_composition_conserves_h_not_u() {
    let source = catalog(200.0);
    let mut gas = original(&source);
    let mass = gas.mass_kg();
    let step = gas.replace_model_with_heat_at_constant_pressure(
        composition(&source, [0.0, 1.0]), 0.0, controls(), || false,
    ).unwrap();
    let target_h_over_r = (h_a(500.0) + 300000.0) * 0.056 / R;
    let target = 2.0 * target_h_over_r /
        (6.5 + (6.5_f64 * 6.5 + 0.004 * target_h_over_r).sqrt());
    close(gas.state().state().temperature, target);
    assert!(target > 500.0 && gas.volume_m3() < 0.01);
    assert!(step.energy.work_by_gas_j < 0.0);
    assert!(step.energy.stored_energy_change_j > 0.0);
    assert!(step.enthalpy_change_j.abs() <= mass * controls().energy_tolerance_j_kg);
    close(step.energy.stored_energy_change_j, -step.energy.work_by_gas_j);
}

#[test]
fn new_composition_needs_only_its_actual_endpoint_not_an_invented_intermediate() {
    let source = catalog(650.0);
    let mut gas = original(&source);
    let heat = gas.mass_kg() * (h_b(900.0) - h_a(500.0));
    let limits = GasEnergySolve { minimum_temperature_k: 650.0, ..controls() };
    gas.replace_model_with_heat_at_constant_pressure(
        composition(&source, [0.0, 1.0]), heat, limits, || false,
    ).unwrap();
    close(gas.state().state().temperature, 900.0);
}

#[test]
fn identical_model_reuses_frozen_heat_bits() {
    let source = catalog(200.0);
    for heat in [0.0, 10.0, -10.0] {
        let mut a = original(&source);
        let mut b = a.clone();
        let model = b.model().clone();
        let expected = a.add_heat_at_constant_pressure(heat, controls(), || false).unwrap();
        let step = b.replace_model_with_heat_at_constant_pressure(model, heat, controls(), || false).unwrap();
        assert_eq!(a, b);
        assert_eq!(step.energy.work_by_gas_j, expected.work_by_gas_j);
        assert_eq!(step.enthalpy_change_j, expected.enthalpy_change_j);
        assert_eq!(step.energy.stored_energy_change_j, expected.stored_energy_change_j);
    }
}

#[test]
fn refused_and_cancelled_isobaric_composition_updates_are_atomic() {
    let source = catalog(200.0);
    let mut gas = original(&source);
    let saved = gas.clone();
    let target = composition(&source, [0.0, 1.0]);
    for (heat, limits) in [(f64::NAN, controls()), (1e12, controls()),
        (0.0, GasEnergySolve { maximum_iterations: 1, ..controls() }),
        (0.0, GasEnergySolve { maximum_temperature_k: 2100.0, ..controls() })]
    {
        assert!(gas.replace_model_with_heat_at_constant_pressure(target.clone(), heat, limits, || false).is_err());
        assert_eq!(gas, saved);
    }
    let mut count = 0;
    gas.clone().replace_model_with_heat_at_constant_pressure(target.clone(), 0.0, controls(),
        || { count += 1; false }).unwrap();
    for stop in [1, count / 2, count] {
        let mut polls = 0;
        assert!(matches!(gas.replace_model_with_heat_at_constant_pressure(target.clone(), 0.0,
            controls(), || { polls += 1; polls == stop }), Err(GasEnergyError::Cancelled)));
        assert_eq!(gas, saved);
    }
}

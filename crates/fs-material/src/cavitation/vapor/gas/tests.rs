use super::*;
use crate::thermochemical_gas::tests::{
    model, model_in_domains, point, species_with_energy_reference, transport,
};
use crate::thermochemical_gas::ThermochemicalGasModel;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    Composition, CompositionBasis, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};

const INITIAL_TEMPERATURE: f64 = 300.0;
const INITIAL_PRESSURE: f64 = 100_000.0;
const INITIAL_VOLUME: f64 = 4.0e-12;

fn controls() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 270.0,
        maximum_temperature_k: 620.0,
        energy_tolerance_j_kg: 1.0e-8,
        maximum_iterations: 64,
    }
}

fn parcel(model: ThermochemicalGasModel) -> ClosedGasParcel {
    ClosedGasParcel::new(
        model,
        point(INITIAL_TEMPERATURE, INITIAL_PRESSURE)
            .with("specimen-lot", 42.0)
            .unwrap(),
        INITIAL_VOLUME,
    )
    .unwrap()
}

fn gas() -> ConservedBubbleGas {
    ConservedBubbleGas::new(parcel(model(vec![0.75, 0.25], 0.001)), controls()).unwrap()
}

fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * actual.abs().max(expected.abs()),
        "actual={actual:e}, expected={expected:e}",
    );
}

#[test]
fn original_nonlinear_capacities_and_fixed_inventory_determine_every_trial() {
    let gas = gas();
    let original = gas.clone();
    let gas_constant = R / (0.75 * 0.028 + 0.25 * 0.004);
    for temperature in [275.0, 300.0, 373.15, 600.0] {
        for volume in [0.5 * INITIAL_VOLUME, INITIAL_VOLUME, 2.0 * INITIAL_VOLUME] {
            let sample = gas.sample(temperature, volume).unwrap();
            let cp_over_r = 0.75 * (3.5 + 0.001 * temperature) + 0.25 * 2.5;
            let expected_cv = (cp_over_r - 1.0) * gas_constant;
            let expected_pressure = INITIAL_PRESSURE * (temperature / INITIAL_TEMPERATURE)
                * (INITIAL_VOLUME / volume);
            assert_eq!(sample.mass_kg, gas.reference_parcel().mass_kg());
            assert_eq!(sample.temperature_k, temperature);
            assert_eq!(sample.volume_m3, volume);
            near(sample.gas_constant_j_kg_k, gas_constant, 2.0e-14);
            near(sample.cv_j_kg_k, expected_cv, 2.0e-14);
            near(sample.cp_j_kg_k, cp_over_r * gas_constant, 2.0e-14);
            near(sample.heat_capacity_j_k, gas.mass_kg() * expected_cv, 2.0e-14);
            near(sample.pressure_pa, expected_pressure, 2.0e-14);
            near(sample.density_kg_m3 * volume, gas.mass_kg(), 2.0e-14);
            near(sample.gamma, cp_over_r / (cp_over_r - 1.0), 2.0e-14);
            near(
                sample.sound_speed_m_s,
                (sample.gamma * gas_constant * temperature).sqrt(),
                2.0e-14,
            );
            let qualified = gas.qualify(temperature, volume).unwrap();
            near(qualified.state().pressure, sample.pressure_pa, 2.0e-14);
            near(qualified.state().density, sample.density_kg_m3, 2.0e-14);
            near(
                qualified.thermodynamics().properties().mass_specific_cv().value(),
                sample.cv_j_kg_k,
                2.0e-14,
            );
        }
    }
    assert_eq!(gas, original);
    assert_eq!(gas.controls(), controls());
}

fn shifted_gas(reference_a: f64, reference_b: f64) -> ConservedBubbleGas {
    let source = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_with_energy_reference("A", 0.028, 3.5, 0.001, reference_a),
                species_with_energy_reference("B", 0.004, 2.5, 0.0, reference_b),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2.0e-5, 0.15)],
        None,
    )
    .unwrap();
    ConservedBubbleGas::new(parcel(source), controls()).unwrap()
}

#[test]
fn caloric_increments_integrate_original_cv_and_ignore_species_formation_datums() {
    let gas = shifted_gas(0.0, 0.0);
    let shifted = shifted_gas(1.0e10, -3.0e9);
    let gas_constant = R / (0.75 * 0.028 + 0.25 * 0.004);
    for (initial, final_temperature) in [(300.0, 400.0), (500.0, 280.0), (300.0, 300.0 + 1.0e-8)] {
        // Independently integrate Cv/R = 2.25 + 0.00075 T, factoring the
        // temperature difference instead of subtracting large stored energies.
        let expected = gas.mass_kg() * gas_constant * (final_temperature - initial)
            * (2.25 + 0.000375 * (final_temperature + initial));
        let change = gas.energy_change_j(initial, final_temperature).unwrap();
        near(change, expected, 3.0e-14);
        assert_eq!(change, shifted.energy_change_j(initial, final_temperature).unwrap());
        near(gas.energy_change_j(final_temperature, initial).unwrap(), -change, 2.0e-14);
        assert_eq!(
            gas.sample(final_temperature, INITIAL_VOLUME).unwrap(),
            shifted.sample(final_temperature, INITIAL_VOLUME).unwrap(),
        );
    }
    assert_eq!(gas.energy_change_j(300.0, 300.0).unwrap(), 0.0);
}

#[test]
fn full_endpoint_qualification_retains_original_sources_and_auxiliary_conditions() {
    let gas = gas();
    let before = gas.clone();
    let temperature = 350.0;
    let volume = 1.25 * INITIAL_VOLUME;
    let sample = gas.sample(temperature, volume).unwrap();
    let qualified = gas.qualify(temperature, volume).unwrap();
    for (current, initial) in qualified
        .components()
        .iter()
        .zip(gas.reference_parcel().state().components())
    {
        assert_eq!(current.card_identity(), initial.card_identity());
        assert_eq!(current.state_coordinate("specimen-lot"), Some(42.0));
        assert_eq!(current.state_coordinate("temperature"), Some(temperature));
        assert_eq!(current.state_coordinate("pressure"), Some(sample.pressure_pa));
        assert_eq!(current.axis_quantities(), gas.reference_parcel().point().axis_quantities());
    }
    assert_eq!(gas, before);
    assert_eq!(gas.reference_parcel().volume_m3(), INITIAL_VOLUME);
    assert_eq!(gas.reference_parcel().state().state().temperature, INITIAL_TEMPERATURE);
}

#[test]
fn source_branch_water_and_transport_refusals_leave_the_inventory_unchanged() {
    let source = model_in_domains(
        vec![0.75, 0.25],
        0.001,
        [200.0, 1000.0],
        [280.0, 380.0],
        [10_000.0, 200_000.0],
    );
    let gas = ConservedBubbleGas::new(parcel(source), controls()).unwrap();
    let before = gas.clone();
    // A numeric trial is not a transport/source acceptance receipt.
    assert!(gas.sample(400.0, INITIAL_VOLUME).is_ok());
    assert!(gas.qualify(400.0, INITIAL_VOLUME).is_err());
    assert!(gas.sample(300.0, 0.1 * INITIAL_VOLUME).is_ok());
    assert!(gas.qualify(300.0, 0.1 * INITIAL_VOLUME).is_err());
    for temperature in [f64::NAN, f64::INFINITY, 269.99, 620.01] {
        assert!(gas.sample(temperature, INITIAL_VOLUME).is_err());
        assert!(gas.energy_change_j(300.0, temperature).is_err());
    }
    for volume in [0.0, -INITIAL_VOLUME, f64::NAN, f64::INFINITY] {
        assert!(gas.sample(300.0, volume).is_err());
    }
    assert_eq!(gas, before);

    let original = model(vec![0.75, 0.25], 0.001);
    let water = SpeciesId::new("A").unwrap();
    let contains_declared_water = ThermochemicalGasModel::new(
        original.thermodynamics().clone(),
        original.transport_materials().to_vec(),
        Some(&water),
    )
    .unwrap();
    assert!(ConservedBubbleGas::new(parcel(contains_declared_water), controls()).is_err());

    let mut outside = controls();
    outside.minimum_temperature_k = 400.0;
    assert!(ConservedBubbleGas::new(parcel(original.clone()), outside).is_err());
    outside.minimum_temperature_k = 270.0;
    outside.maximum_temperature_k = 2500.0;
    assert!(ConservedBubbleGas::new(parcel(original), outside).is_err());
}

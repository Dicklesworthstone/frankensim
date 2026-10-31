use super::*;
use crate::phase::condensation::{WaterCondensationModel, WetGasParcel};
use crate::thermochemical_gas::energy::IsochoricGasParcel;
use crate::thermochemical_gas::tests::{point, species_with_energy_reference, transport};
use fs_qty::Temperature;
use fs_thermochem::{SpeciesId, mixture::FrozenIdealGasMixtureModelV1, water::SaturatedWaterState};

fn source() -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_with_energy_reference("carrier", 0.028, 3.5, 0.0, 0.0),
                species_with_energy_reference("water", 0.01801528, 4.0, 0.001, -13.4e6),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![0.9, 0.1]).unwrap(),
        )
        .unwrap(),
        vec![
            transport("carrier", 1.7e-5, 0.025),
            transport("water", 9e-6, 0.019),
        ],
        Some(&SpeciesId::new("water").unwrap()),
    )
    .unwrap()
}

#[test]
fn g0_iapws_condensation_endpoint_feeds_the_original_shared_field_without_a_second_phase_curve() {
    let source = source();
    let gas = IsochoricGasParcel::new(source.clone(), point(300.0, 101325.0), 0.001).unwrap();
    let law = WaterCondensationModel::new(source, SpeciesId::new("water").unwrap(), 300.0).unwrap();
    let wet = WetGasParcel::new(law, gas, 0.0)
        .unwrap()
        .at_temperature(300.0, || false)
        .unwrap()
        .parcel;
    assert!(wet.liquid_mass_kg() > 0.0);
    let gas = wet.gas();
    let primitive = PrimitiveGasState {
        temperature_k: wet.temperature_k(),
        composition: gas.model().thermodynamics().mole_composition().clone(),
    };
    let field = AxialGasField::new(
        vec![0.0, 0.5],
        gas.state().state().pressure,
        vec![primitive.clone(), primitive],
        gas.model().clone(),
        QueryPoint::new(),
        ThermochemicalPhase::LiquidWaterIapwsLimit,
    )
    .unwrap();
    for x in [0.0, 0.1, 0.25, 0.5] {
        assert_eq!(field.state_at(x).unwrap(), *gas.state());
    }
    assert_eq!(
        field.sample_states(&[0.0, 0.5]).unwrap(),
        vec![*gas.state().state(); 2]
    );
    assert!(
        gas.model()
            .evaluate_with_liquid_water_limit(gas.point())
            .is_err(),
        "the legacy Buck curve is distinct and must not silently replace IAPWS"
    );
}

#[test]
fn g0_iapws_admission_refuses_real_excess_water_and_unsupported_condensed_phase() {
    let source = source();
    for (t, p) in [(300.0, 101325.0), (273.15, 101325.0), (300.0, 300000.0)] {
        assert!(
            source
                .evaluate_with_iapws_liquid_water_limit(&point(t, p))
                .is_err()
        );
    }
    let temperature = 300.0;
    let pressure = 101325.0;
    let sat = SaturatedWaterState::at_temperature(Temperature::new(temperature)).unwrap();
    let fraction = sat.pressure().value() / pressure;
    let at = |water| {
        source
            .with_composition(
                Composition::new(CompositionBasis::MoleFraction, vec![1.0 - water, water]).unwrap(),
            )
            .unwrap()
    };
    assert!(
        at(fraction)
            .evaluate_with_iapws_liquid_water_limit(&point(temperature, pressure))
            .is_ok()
    );
    assert!(
        at(fraction * (1.0 + 1e-10))
            .evaluate_with_iapws_liquid_water_limit(&point(temperature, pressure))
            .is_err()
    );
}

#[test]
fn g0_a_pair_of_saturated_stations_does_not_silently_clamp_a_supersaturated_interpolant() {
    let pressure = 101325.0;
    let states = [280.0, 320.0].map(|temperature_k| {
        let water = SaturatedWaterState::at_temperature(Temperature::new(temperature_k))
            .unwrap()
            .pressure()
            .value()
            / pressure;
        PrimitiveGasState {
            temperature_k,
            composition: Composition::new(CompositionBasis::MoleFraction, vec![1.0 - water, water])
                .unwrap(),
        }
    });
    let field = AxialGasField::new(
        vec![0.0, 1.0],
        pressure,
        states.to_vec(),
        source(),
        QueryPoint::new(),
        ThermochemicalPhase::LiquidWaterIapwsLimit,
    )
    .unwrap();
    assert!(field.state_at(0.0).is_ok());
    assert!(field.state_at(1.0).is_ok());
    assert!(
        field
            .state_at(0.5)
            .unwrap_err()
            .to_string()
            .contains("supersaturated")
    );
}

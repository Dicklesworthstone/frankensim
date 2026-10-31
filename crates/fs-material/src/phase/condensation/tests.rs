use super::*;
use crate::thermochemical_gas::tests::{point, species_with_energy_reference, transport};
use fs_thermochem::{Composition, CompositionBasis, mixture::FrozenIdealGasMixtureModelV1};

fn source(water_fraction: f64, reference_shift: f64) -> ThermochemicalGasModel {
    let species = vec![
        species_with_energy_reference("carrier", 0.028, 3.5, 0.0, 0.0),
        species_with_energy_reference("water", 0.01801528, 4.0, 0.001, -13.4e6 + reference_shift),
    ];
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            species,
            Composition::new(
                CompositionBasis::MoleFraction,
                vec![1.0 - water_fraction, water_fraction],
            )
            .unwrap(),
        )
        .unwrap(),
        vec![
            transport("carrier", 1.7e-5, 0.025),
            transport("water", 9.0e-6, 0.019),
        ],
        Some(&SpeciesId::new("water").unwrap()),
    )
    .unwrap()
}
fn wet(t: f64, fraction: f64, liquid: f64, shift: f64) -> WetGasParcel {
    let source = source(fraction, shift);
    let model =
        WaterCondensationModel::new(source.clone(), SpeciesId::new("water").unwrap(), t).unwrap();
    let gas = IsochoricGasParcel::new(
        source,
        point(t, 101325.0).with("specimen-lot", 17.0).unwrap(),
        0.001,
    )
    .unwrap();
    WetGasParcel::new(model, gas, liquid).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn g0_cooling_condenses_against_partial_pressure_and_retains_every_carrier() {
    let original = wet(310.0, 0.06, 0.0, 0.0);
    let before = original.clone();
    let step = original.at_temperature(290.0, || false).unwrap();
    let next = &step.parcel;
    assert!(step.condensed_mass_kg > 0.0);
    assert_eq!(next.total_water_mass_kg(), original.total_water_mass_kg());
    let carrier = SpeciesId::new("carrier").unwrap();
    assert_eq!(
        next.species_masses_kg()[&carrier],
        original.species_masses_kg()[&carrier]
    );
    let water_r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.01801528;
    let saturation = SaturatedWaterState::at_temperature(Temperature::new(290.0)).unwrap();
    let expected = saturation.pressure().value() * 0.001 / (water_r * 290.0);
    near(next.vapor_mass_kg(), expected, expected * 4.0e-15);
    near(
        next.vapor_mass_kg() + next.liquid_mass_kg(),
        original.total_water_mass_kg(),
        1e-18,
    );
    near(step.mass_residual_kg, 0.0, original.gas().mass_kg() * 1e-13);
    near(step.saturation_residual_pa, 0.0, 1e-9);
    assert_eq!(next.gas().volume_m3(), original.gas().volume_m3());
    assert_eq!(next.gas().point().axes().get("specimen-lot"), Some(&17.0));
    assert!(next.gas().state().state().pressure < original.gas().state().state().pressure);
    assert_ne!(
        next.gas().state().state().sound_speed,
        original.gas().state().state().sound_speed
    );
    assert_eq!(
        original, before,
        "proposing an endpoint must not mutate the accepted gas"
    );
}

#[test]
fn g0_isothermal_condensation_releases_internal_energy_without_double_counting_latent_heat() {
    let original = wet(300.0, 0.10, 0.0, 0.0);
    let step = original.at_temperature(300.0, || false).unwrap();
    let sat = original.liquid_saturation();
    let latent_h = sat.latent_enthalpy().value();
    let latent_u = latent_h - UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.01801528 * 300.0
        + sat.pressure().value() / sat.liquid_density().value();
    near(
        step.latent_enthalpy_j,
        step.condensed_mass_kg * latent_h,
        1e-10,
    );
    near(
        step.total_internal_energy_change_j,
        -step.condensed_mass_kg * latent_u,
        1e-10,
    );
    assert!(step.total_internal_energy_change_j < 0.0);
    assert!(
        (step.total_internal_energy_change_j + step.latent_enthalpy_j).abs() > 1.0,
        "rigid gas storage uses latent internal energy, not latent enthalpy"
    );
    near(
        step.gas_internal_energy_change_j + step.liquid_internal_energy_change_j,
        step.total_internal_energy_change_j,
        1e-9,
    );
    near(
        step.parcel.total_internal_energy_j() - original.total_internal_energy_j(),
        step.total_internal_energy_change_j,
        1e-9,
    );
}

#[test]
fn g0_evaporation_exhausts_only_the_finite_film_then_remains_unsaturated() {
    let original = wet(300.0, 0.0001, 1e-6, 0.0);
    let step = original.at_temperature(310.0, || false).unwrap();
    assert_eq!(step.parcel.liquid_mass_kg(), 0.0);
    assert_eq!(step.parcel.vapor_mass_kg(), original.total_water_mass_kg());
    near(step.condensed_mass_kg, -1e-6, 1e-20);
    assert!(step.saturation_residual_pa < 0.0);
    assert!(step.total_internal_energy_change_j > 0.0);
    let repeat = step.parcel.at_temperature(310.0, || false).unwrap();
    assert_eq!(repeat.condensed_mass_kg, 0.0);
    assert_eq!(repeat.total_internal_energy_change_j, 0.0);
}

#[test]
fn g0_exactly_dry_active_gas_can_recover_water_from_the_retained_catalog() {
    let catalog = source(0.05, 0.0);
    let model =
        WaterCondensationModel::new(catalog.clone(), SpeciesId::new("water").unwrap(), 300.0)
            .unwrap();
    let dry = catalog
        .with_composition(Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0]).unwrap())
        .unwrap();
    assert!(dry.water_species().is_none());
    let gas = IsochoricGasParcel::new(dry, point(300.0, 101325.0), 0.001).unwrap();
    let wet = WetGasParcel::new(model, gas, 1e-4).unwrap();
    let step = wet.at_temperature(300.0, || false).unwrap();
    assert!(step.condensed_mass_kg < 0.0);
    assert!(step.parcel.liquid_mass_kg() > 0.0);
    assert_eq!(
        step.parcel.gas().model().water_species(),
        catalog.water_species()
    );
    assert_eq!(
        step.parcel.gas().model().thermodynamics().components(),
        catalog.thermodynamics().components()
    );
    near(step.saturation_residual_pa, 0.0, 1e-9);
}

#[test]
fn g3_phase_heat_is_invariant_to_large_chemical_formation_datum_shifts() {
    let original = wet(310.0, 0.06, 2e-5, 0.0);
    let shifted = wet(310.0, 0.06, 2e-5, 1e15);
    let a = original.at_temperature(290.0, || false).unwrap();
    let b = shifted.at_temperature(290.0, || false).unwrap();
    assert_eq!(a.condensed_mass_kg, b.condensed_mass_kg);
    near(
        a.total_internal_energy_change_j,
        b.total_internal_energy_change_j,
        1e-10,
    );
    assert_ne!(
        a.gas_internal_energy_change_j,
        b.gas_internal_energy_change_j
    );
    assert_ne!(
        a.liquid_internal_energy_change_j,
        b.liquid_internal_energy_change_j
    );
}

#[test]
fn g5_repeated_cooling_and_evaporation_retain_canonical_species_without_drift() {
    let original = wet(310.0, 0.04, 2e-5, 0.0);
    let mut parcel = original.clone();
    let carrier = SpeciesId::new("carrier").unwrap();
    for t in [285.0, 310.0, 285.0, 310.0, 285.0, 310.0] {
        let a = parcel.at_temperature(t, || false).unwrap();
        let b = parcel.at_temperature(t, || false).unwrap();
        assert_eq!(a, b);
        parcel = a.parcel;
        assert_eq!(parcel.total_water_mass_kg(), original.total_water_mass_kg());
        assert_eq!(
            parcel.species_masses_kg()[&carrier],
            original.species_masses_kg()[&carrier]
        );
    }
}

#[test]
fn g4_unrepresentable_phase_inventory_and_transfer_are_refused_atomically() {
    let original = wet(300.0, 0.10, 0.0, 0.0);
    assert!(WetGasParcel::new(original.model().clone(), original.gas().clone(), 1e13).is_err());
    assert!(WetGasParcel::new(original.model().clone(), original.gas().clone(), 1e-30).is_err());

    // An exactly dry gas has no vapor to lose in the initial total, but a film
    // whose next evaporation increment cannot be represented must still refuse.
    let catalog = source(0.05, 0.0);
    let model =
        WaterCondensationModel::new(catalog.clone(), SpeciesId::new("water").unwrap(), 300.0)
            .unwrap();
    let dry = catalog
        .with_composition(Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0]).unwrap())
        .unwrap();
    let gas = IsochoricGasParcel::new(dry, point(300.0, 101325.0), 0.001).unwrap();
    let original = WetGasParcel::new(model, gas, 1e13).unwrap();
    let before = original.clone();
    assert!(matches!(
        original.at_temperature(300.0, || false),
        Err(WaterCondensationError::Invalid(_))
    ));
    assert_eq!(original, before);
}

#[test]
fn g4_cancellation_source_mismatch_and_unsupported_phase_domains_are_atomic() {
    let original = wet(310.0, 0.06, 0.0, 0.0);
    let before = original.clone();
    assert_eq!(
        original.at_temperature(290.0, || true),
        Err(WaterCondensationError::Cancelled)
    );
    let mut checks = 0;
    assert_eq!(
        original.at_temperature(290.0, || {
            checks += 1;
            checks >= 3
        }),
        Err(WaterCondensationError::Cancelled)
    );
    for t in [273.15, 323.16, f64::NAN, f64::INFINITY] {
        assert!(original.at_temperature(t, || false).is_err());
    }
    let foreign =
        IsochoricGasParcel::new(source(0.06, 1e4), point(310.0, 101325.0), 0.001).unwrap();
    assert!(WetGasParcel::new(original.model().clone(), foreign, 0.0).is_err());
    let overpressure =
        IsochoricGasParcel::new(source(0.06, 0.0), point(310.0, 300000.0), 0.001).unwrap();
    assert!(WetGasParcel::new(original.model().clone(), overpressure, 0.0).is_err());
    let wrong_water = SpeciesId::new("carrier").unwrap();
    let template = source(0.06, 0.0);
    let wrong = ThermochemicalGasModel::new(
        template.thermodynamics().clone(),
        template.transport_materials().to_vec(),
        Some(&wrong_water),
    )
    .unwrap();
    assert!(WaterCondensationModel::new(wrong, wrong_water, 300.0).is_err());
    assert_eq!(original, before);
}

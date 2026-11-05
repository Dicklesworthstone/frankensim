//! G0/G1 source ownership and independent ideal-plasma scale checks.

use super::*;
use crate::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, InitialStatePolicy,
    MolarEnergyQuantityV1, MolarThermalQuantityV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1,
    NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReferenceEquationOfStateV1, StandardStateConventionV1, StandardStatePhaseV1,
    required_nasa9_parameter_dims,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, LawId, LawParameter, Provenance};
use fs_qty::{MolarMass, Pressure};
use std::collections::BTreeMap;

const NAMES: [&str; 4] = ["ion_d", "ion_t", "ash_he", "electron"];

fn id(name: &str) -> SpeciesId { SpeciesId::new(name).unwrap() }

fn source(name: &str, mass: f64) -> Nasa9StandardStateModelV1 {
    let coefficients = [0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let mut parameters = BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(coefficients) {
        parameters.insert(name.into(), LawParameter { value, dims });
    }
    parameters.insert("reference_pressure".into(), LawParameter { value: 1e5, dims: Pressure::DIMS });
    let card = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()),
        law_version: NASA9_LAW_VERSION_V1,
        parameters,
        state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 1e4, 2e9),
        sources: vec![],
        provenance: Provenance {
            source: "manufactured ideal monatomic sensible calorics; no measured plasma claim".into(),
            license: "analytic test fixture".into(),
            artifact: None,
        },
    };
    Nasa9StandardStateModelV1::new(
        id(name), MolarMass::new(mass * AVOGADRO_PER_MOL),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas, ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(1e5), ElementalReferenceIdV1::new("fusion-sensible-test").unwrap(),
        ).unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    ).unwrap()
}

fn sources() -> Vec<Nasa9StandardStateModelV1> {
    [NuclearSpeciesV1::Deuterium.mass_kg(), NuclearSpeciesV1::Tritium.mass_kg(),
        NuclearSpeciesV1::Helium4.mass_kg(), ELECTRON_MASS_KG]
        .into_iter().zip(NAMES).map(|(mass, name)| source(name, mass)).collect()
}

fn bank(models: Vec<Nasa9StandardStateModelV1>, charges: Vec<i128>) -> IdealGasEquilibriumV1 {
    let names: Vec<_> = NAMES.map(id).into();
    IdealGasEquilibriumV1::new(
        ElementalMatrix::new(vec![ElementId::new("H").unwrap(), ElementId::new("He").unwrap()],
            names.clone(), vec![vec![1, 1, 0, 0], vec![0, 0, 1, 0]]).unwrap(),
        ChargeVector::new(names, charges).unwrap(), models,
    ).unwrap()
}

fn bind(bank: IdealGasEquilibriumV1) -> Result<DtFusionPlasmaV1, DtFusionPlasmaErrorV1> {
    DtFusionPlasmaV1::new(bank, id(NAMES[0]), id(NAMES[1]), id(NAMES[2]))
}

fn fixture() -> DtFusionPlasmaV1 { bind(bank(sources(), vec![1, 1, 2, -1])).unwrap() }

fn amounts(plasma: &DtFusionPlasmaV1, heavy: [f64; 3]) -> [Amount; 4] {
    let mut result = [Amount::new(0.0); 4];
    for (index, amount) in plasma.indices().into_iter().zip(heavy) {
        result[index] = Amount::new(amount);
    }
    result[plasma.electron_index()] = Amount::new(heavy[0] + heavy[1] + 2.0 * heavy[2]);
    result
}

fn state(plasma: &DtFusionPlasmaV1, heavy: [f64; 3], th: f64, te: f64) -> FrozenIonizedGasStateV1 {
    FrozenIonizedGasStateV1::new(plasma.thermodynamics(), &amounts(plasma, heavy), Volume::new(2.0),
        Temperature::new(th), Temperature::new(te), 1e-12, || false).unwrap()
}

fn limits() -> DtFusionPlasmaLimitsV1 {
    DtFusionPlasmaLimitsV1 { maximum_coulomb_coupling: 0.01,
        maximum_thermal_phase_space_density: 0.001, maximum_thermal_energy_rest_energy_ratio: 0.05 }
}

fn changed_source(model: &Nasa9StandardStateModelV1, mut edit: impl FnMut(&mut ConstitutiveModelCard)) -> Nasa9StandardStateModelV1 {
    let mut card = model.regions()[0].card().clone();
    edit(&mut card);
    Nasa9StandardStateModelV1::new(model.species().clone(), model.molar_mass(), model.convention().clone(),
        vec![Nasa9RegionV1::from_card(card).unwrap()]).unwrap()
}

#[test]
fn g0_explicit_isotope_mapping_retains_original_source_bank_and_state() {
    let original = bank(sources(), vec![1, 1, 2, -1]);
    let plasma = bind(original.clone()).unwrap();
    assert_eq!(plasma.indices(), [2, 3, 0]);
    assert_eq!(plasma.electron_index(), 1);
    assert!(core::ptr::eq(original.models().as_ptr(), plasma.thermodynamics().models().as_ptr()));
    let point = state(&plasma, [0.001, 0.002, 0.0005], 1e8, 2e7);
    let before = point.clone();
    plasma.admit_state(&point, || false).unwrap();
    plasma.admit_state_with_limits(&point, limits(), || false).unwrap();
    assert_eq!(point, before);
    let r = crate::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
    assert!((point.pressure().value() / (r * (0.0035 * 1e8 + 0.004 * 2e7) / 2.0) - 1.0).abs() < 1e-15);
    assert!((point.heavy_internal_energy().value() / (1.5 * r * 0.0035 * 1e8) - 1.0).abs() < 1e-15);
}

#[test]
fn g0_rejects_isotope_aliases_swapped_masses_and_atomic_mass_sources() {
    let original = bank(sources(), vec![1, 1, 2, -1]);
    for names in [[NAMES[1], NAMES[0], NAMES[2]], [NAMES[0], NAMES[0], NAMES[2]], ["missing", NAMES[1], NAMES[2]]] {
        assert!(DtFusionPlasmaV1::new(original.clone(), id(names[0]), id(names[1]), id(names[2])).is_err());
    }
    let mut models = sources();
    models[0] = source(NAMES[0], NuclearSpeciesV1::Deuterium.mass_kg() + ELECTRON_MASS_KG);
    assert!(bind(bank(models, vec![1, 1, 2, -1])).is_err());
    assert!(bind(bank(sources(), vec![1, 1, 1, -1])).is_err());
}

#[test]
fn g0_nuclear_energy_cannot_be_duplicated_in_original_calorics() {
    for (parameter, value) in [("a7", 1.0), ("a2", 3.5), ("a3", 1e-8)] {
        let mut models = sources();
        models[2] = changed_source(&models[2], |card| card.parameters.get_mut(parameter).unwrap().value = value);
        assert!(bind(bank(models, vec![1, 1, 2, -1])).is_err());
    }
    let mut models = sources();
    models[0] = models[0].with_caloric_continuity(Temperature::new(1e8),
        MolarEnergyQuantityV1::new(0.0), MolarThermalQuantityV1::new(0.0)).unwrap();
    assert!(bind(bank(models, vec![1, 1, 2, -1])).is_err());
}

#[test]
fn g0_unchanged_state_values_do_not_hide_changed_source_identity() {
    let plasma = fixture();
    let original = state(&plasma, [0.001, 0.002, 0.0], 1e8, 2e7);
    let mut models = sources();
    models[0] = changed_source(&models[0], |card| card.provenance.source.push_str(" alternate source"));
    let other = bind(bank(models, vec![1, 1, 2, -1])).unwrap();
    assert!(other.admit_state(&original, || false).is_err());
    assert_eq!(plasma.admit_state(&original, || true), Err(DtFusionPlasmaErrorV1::Cancelled));
}

#[test]
fn g0_possible_ash_needs_original_source_support_even_when_absent() {
    let mut models = sources();
    models[2] = changed_source(&models[2], |card| card.validity = ValidityDomain::unconstrained().with("T", 1e4, 1e7));
    let plasma = bind(bank(models, vec![1, 1, 2, -1])).unwrap();
    let initial = state(&plasma, [0.001, 0.002, 0.0], 1e8, 2e7);
    assert!(matches!(plasma.admit_state(&initial, || false), Err(DtFusionPlasmaErrorV1::Source { .. })));
    let initial = state(&plasma, [0.001, 0.002, 0.0], 1e6, 2e7);
    plasma.admit_state(&initial, || false).unwrap();
}

#[test]
fn g1_physical_scales_match_independent_direct_dimensional_formulas() {
    let plasma = fixture();
    let initial = state(&plasma, [0.001, 0.002, 0.0005], 1e8, 2e7);
    let actual = plasma.admit_state_with_limits(&initial, limits(), || false).unwrap();
    let densities: Vec<_> = (0..4).map(|i| initial.number_density_per_m3(i).unwrap()).collect();
    let total_density: f64 = densities.iter().sum();
    let pi = core::f64::consts::PI;
    let radius = (3.0 / (4.0 * pi * total_density)).cbrt();
    let mut expected = [0.0_f64; 3];
    for (i, n) in densities.into_iter().enumerate() {
        let mass = plasma.thermodynamics().models()[i].molar_mass().value() / AVOGADRO_PER_MOL;
        let temperature = initial.species_temperature(i).unwrap().value();
        let z = plasma.thermodynamics().charge_vector().get(i).unwrap() as f64;
        expected[0] = expected[0].max(z * z * ELEMENTARY_CHARGE_C.powi(2)
            / (4.0 * pi * VACUUM_PERMITTIVITY_F_M * radius * BOLTZMANN_J_PER_K * temperature));
        expected[1] = expected[1].max(n * (PLANCK_J_S / (2.0 * pi * mass * BOLTZMANN_J_PER_K * temperature).sqrt()).powi(3));
        expected[2] = expected[2].max(BOLTZMANN_J_PER_K * temperature / (mass * SPEED_OF_LIGHT_M_S.powi(2)));
    }
    for (got, want) in [actual.maximum_coulomb_coupling, actual.maximum_thermal_phase_space_density,
        actual.maximum_thermal_energy_rest_energy_ratio].into_iter().zip(expected) {
        assert!((got / want - 1.0).abs() < 1e-12, "actual={got}, expected={want}");
    }
    let mut restrictive = limits();
    restrictive.maximum_coulomb_coupling = 0.5 * expected[0];
    assert_eq!(plasma.admit_state_with_limits(&initial, restrictive, || false), Err(DtFusionPlasmaErrorV1::Domain("Coulomb coupling")));
    restrictive = limits();
    restrictive.maximum_thermal_phase_space_density = 0.5 * expected[1];
    assert_eq!(plasma.admit_state_with_limits(&initial, restrictive, || false), Err(DtFusionPlasmaErrorV1::Domain("thermal phase-space density")));
    restrictive = limits();
    restrictive.maximum_thermal_energy_rest_energy_ratio = 0.5 * expected[2];
    assert_eq!(plasma.admit_state_with_limits(&initial, restrictive, || false), Err(DtFusionPlasmaErrorV1::Domain("nonrelativistic thermal energy")));
}

#[test]
fn g0_trial_admission_refuses_invalid_inventory_controls_and_cancellation() {
    let plasma = fixture();
    let original = amounts(&plasma, [0.001, 0.002, 0.0]);
    let admit = |amounts: &[Amount], volume, physical_limits| plasma.admit_point(Temperature::new(1e8),
        Temperature::new(2e7), amounts, Volume::new(volume), physical_limits, || false);
    admit(&original, 2.0, limits()).unwrap();
    for bad in [-1.0, f64::INFINITY, f64::NAN] {
        let mut changed = original;
        changed[plasma.indices()[0]] = Amount::new(bad);
        assert!(admit(&changed, 2.0, limits()).is_err());
    }
    assert!(admit(&[Amount::new(0.0); 4], 2.0, limits()).is_err());
    assert!(admit(&original[..3], 2.0, limits()).is_err());
    assert!(admit(&original, 0.0, limits()).is_err());
    let mut changed = original;
    changed[plasma.electron_index()] = Amount::new(0.0);
    assert!(admit(&changed, 2.0, limits()).is_err());
    let mut invalid_limits = limits();
    invalid_limits.maximum_coulomb_coupling = f64::NAN;
    assert!(admit(&original, 2.0, invalid_limits).is_err());
    assert_eq!(plasma.admit_point(Temperature::new(1e8), Temperature::new(2e7), &original,
        Volume::new(2.0), limits(), || true), Err(DtFusionPlasmaErrorV1::Cancelled));
    // Fuel exhaustion remains an admitted source state; no artificial seeds.
    admit(&amounts(&plasma, [0.0, 0.0, 0.001]), 2.0, limits()).unwrap();
}

#[test]
fn g0_contiguous_original_regions_are_admitted_but_gaps_and_entropy_jumps_are_not() {
    for (second_start, entropy_jump) in [(1e7, false), (2e7, false), (1e7, true)] {
        let mut models = sources();
        let mut lower = models[2].regions()[0].card().clone();
        lower.validity = ValidityDomain::unconstrained().with("T", 1e4, 1e7);
        let mut upper = lower.clone();
        upper.validity = ValidityDomain::unconstrained().with("T", second_start, 2e9);
        if entropy_jump { upper.parameters.get_mut("a8").unwrap().value = 1.0; }
        models[2] = Nasa9StandardStateModelV1::new(models[2].species().clone(), models[2].molar_mass(),
            models[2].convention().clone(), vec![Nasa9RegionV1::from_card(lower).unwrap(),
                Nasa9RegionV1::from_card(upper).unwrap()]).unwrap();
        let result = bind(bank(models, vec![1, 1, 2, -1]));
        if entropy_jump { assert!(result.is_err()); continue; }
        let plasma = result.unwrap();
        let result = plasma.admit_temperature_interval(Temperature::new(1e6), Temperature::new(1e8),
            Temperature::new(1e5), Temperature::new(1e7), || false);
        assert_eq!(result.is_ok(), second_start == 1e7);
    }
}

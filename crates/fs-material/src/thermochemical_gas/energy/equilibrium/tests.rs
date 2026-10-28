//! Manufactured isomers test coupling and conservation, not propellant accuracy.
mod packs;

use super::*;
use crate::thermochemical_gas::SpeciesTransportMaterial;
use crate::thermochemical_gas::energy::GasEnergySolve;
use crate::thermochemical_gas::energy::parcel::transfer::nozzle::rocket::{
    RocketNozzle, RocketNozzleSolve,
};
use crate::thermochemical_gas::tests::{
    point, species, species_with_energy_reference, transport, transport_in_domain,
};
use fs_thermochem::equilibrium::composition::EquilibriumSolveV1;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, Nasa9RegionV1, Nasa9StandardStateModelV1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};

fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}

fn controls() -> GasEquilibriumChamberSolve {
    GasEquilibriumChamberSolve {
        equilibrium: EquilibriumEnthalpySolveV1 {
            equilibrium: EquilibriumSolveV1::default(),
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1500.0,
            enthalpy_tolerance_j: 1e-6,
            maximum_temperature_iterations: 128,
        },
        mass_tolerance_kg: 1e-12,
    }
}

pub(super) fn fixture(
    models: Vec<Nasa9StandardStateModelV1>,
    transport: Vec<SpeciesTransportMaterial>,
    elements: &[&str],
    rows: Vec<Vec<i128>>,
) -> (ThermochemicalGasModel, IdealGasEquilibriumV1) {
    let ids: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    let n = models.len();
    let mut fractions = vec![1.0 / n as f64; n];
    close_generated_simplex(&mut fractions).unwrap();
    let bank = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            models.clone(),
            Composition::new(CompositionBasis::MoleFraction, fractions).unwrap(),
        )
        .unwrap(),
        transport,
        None,
    )
    .unwrap();
    let elemental = ElementalMatrix::new(
        elements
            .iter()
            .map(|e| ElementId::new(*e).unwrap())
            .collect(),
        ids.clone(),
        rows,
    )
    .unwrap();
    let chemistry = IdealGasEquilibriumV1::new(
        elemental,
        ChargeVector::new(ids, vec![0; n]).unwrap(),
        models,
    )
    .unwrap();
    (bank, chemistry)
}

pub(super) fn isomers() -> (ThermochemicalGasModel, IdealGasEquilibriumV1) {
    fixture(
        vec![
            species("A", 0.01, 4.0, 0.0),
            species_with_energy_reference("B", 0.01, 4.0, 0.0, -2.4e6),
        ],
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2.0e-5, 0.15)],
        &["H"],
        vec![vec![1, 1]],
    )
}

fn feed() -> Vec<(SpeciesId, Amount)> {
    vec![(id("A"), Amount::new(2.0))]
}

#[test]
fn exothermic_chamber_preserves_sources_and_drives_the_existing_frozen_nozzle() {
    let (bank, chemistry) = isomers();
    let old = bank.clone();
    let context = point(300.0, 3e6).with("sample_id", 7.0).unwrap();
    let result = bank
        .equilibrium_chamber_hp(
            &chemistry,
            &context,
            &feed(),
            Temperature::new(300.0),
            controls(),
            || false,
        )
        .unwrap();
    let gas = result.parcel();
    let t = gas.state().state().temperature;
    // Independent equilibrium fraction and energy root for A <-> B with
    // equal constant cp and a -24 kJ/mol product formation enthalpy.
    let (mut lo, mut hi) = (300.0, 1500.0);
    for _ in 0..70 {
        let trial = 0.5 * (lo + hi);
        let y_b = 1.0 / (1.0 + (-24000.0 / (R * trial)).exp());
        if 4.0 * R * (trial - 300.0) - 24000.0 * y_b < 0.0 {
            lo = trial;
        } else {
            hi = trial;
        }
    }
    assert!((t - 0.5 * (lo + hi)).abs() < 1e-7);
    assert!((gas.mass_kg() - 0.02).abs() < controls().mass_tolerance_kg);
    assert_eq!(
        gas.volume_m3(),
        result.equilibrium().mass().value() / gas.state().state().density
    );
    assert_eq!(gas.point().axes()["sample_id"], 7.0);
    assert_eq!(gas.state().state().pressure, 3e6);
    assert_eq!(
        gas.model().thermodynamics().components(),
        chemistry.models()
    );
    assert_eq!(
        gas.model().transport_materials(),
        bank.transport_materials()
    );
    assert_eq!(
        gas.model()
            .thermodynamics()
            .mole_composition()
            .fractions()
            .iter()
            .sum::<f64>(),
        1.0
    );
    assert!(
        result.balance().enthalpy_residual_j.abs() < controls().equilibrium.enthalpy_tolerance_j
    );
    assert!(
        result.balance().freeze_enthalpy_residual_j.abs()
            < controls().equilibrium.enthalpy_tolerance_j
    );
    assert!(result.balance().maximum_mole_fraction_correction <= 4.0 * f64::EPSILON);
    assert!((result.balance().parcel_enthalpy_j - 2.0 * 4.0 * R * 300.0).abs() < 1e-6);
    let flow = gas
        .rocket_nozzle_flow(
            0.0,
            RocketNozzle {
                throat_area_m2: 1e-4,
                exit_area_m2: 2e-4,
            },
            RocketNozzleSolve {
                nozzle: GasEnergySolve {
                    minimum_temperature_k: 200.0,
                    maximum_temperature_k: 1500.0,
                    energy_tolerance_j_kg: 1e-6,
                    maximum_iterations: 128,
                },
                relative_mass_flow_tolerance: 1e-10,
            },
            || false,
        )
        .unwrap();
    assert!(flow.exit_mach_number > 1.0 && flow.thrust_n > 0.0);
    assert!(flow.exit().state().temperature < t);
    assert_eq!(
        flow.exit().thermodynamics().properties().molar_mass(),
        gas.state().thermodynamics().properties().molar_mass()
    );
    assert!(flow.energy_residual_j_kg.abs() < 2e-6);
    assert!(flow.relative_mass_flow_residual.abs() < 1e-10);
    assert_eq!(bank, old);
}

#[test]
fn separately_heated_streams_close_one_total_enthalpy_without_reaction_heat_duplication() {
    let (bank, chemistry) = isomers();
    let first = [(id("A"), Amount::new(1.5))];
    let second = [(id("B"), Amount::new(0.5))];
    let target = chemistry
        .feed_enthalpy(Temperature::new(300.0), &first, || false)
        .unwrap()
        + chemistry
            .feed_enthalpy(Temperature::new(700.0), &second, || false)
            .unwrap();
    let mut feed = vec![first[0].clone(), second[0].clone()];
    let result = bank
        .equilibrium_chamber_hp_enthalpy(
            &chemistry,
            &point(250.0, 2e6),
            &feed,
            target,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(result.equilibrium().target_enthalpy(), Some(target));
    assert!((result.balance().parcel_enthalpy_j - target.value()).abs() < 1e-6);
    assert!(result.balance().mass_residual_kg.abs() < 1e-12);
    // Catalog fractions and input ordering do not select chamber chemistry.
    let other = bank
        .with_composition(Composition::new(CompositionBasis::MoleFraction, vec![0.9, 0.1]).unwrap())
        .unwrap();
    feed.reverse();
    let repeated = other
        .equilibrium_chamber_hp_enthalpy(
            &chemistry,
            &point(250.0, 2e6),
            &feed,
            target,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(result, repeated);
}

#[test]
fn exact_source_identity_and_complete_species_bank_are_required() {
    let (bank, chemistry) = isomers();
    let changed = fixture(
        vec![species("A", 0.01, 4.01, 0.0), chemistry.models()[1].clone()],
        bank.transport_materials().to_vec(),
        &["H"],
        vec![vec![1, 1]],
    )
    .0;
    assert_eq!(
        changed.equilibrium_chamber_hp(
            &chemistry,
            &point(300.0, 1e5),
            &feed(),
            Temperature::new(300.0),
            controls(),
            || false,
        ),
        Err(GasEquilibriumChamberError::SourceMismatch(id("A")))
    );
    // Same values with a different provenance card also must not be substituted.
    let model = &chemistry.models()[0];
    let mut card = model.regions()[0].card().clone();
    card.provenance.source = "different authored source".into();
    let different_source = Nasa9StandardStateModelV1::new(
        model.species().clone(),
        model.molar_mass(),
        model.convention().clone(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap();
    let changed = fixture(
        vec![different_source, chemistry.models()[1].clone()],
        bank.transport_materials().to_vec(),
        &["H"],
        vec![vec![1, 1]],
    )
    .0;
    assert!(matches!(
        changed.equilibrium_chamber_hp(
            &chemistry,
            &point(300.0, 1e5),
            &feed(),
            Temperature::new(300.0),
            controls(),
            || false,
        ),
        Err(GasEquilibriumChamberError::SourceMismatch(_))
    ));
    let missing = bank
        .with_composition(Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0]).unwrap())
        .unwrap();
    assert!(matches!(
        missing.equilibrium_chamber_hp(
            &chemistry,
            &point(300.0, 1e5),
            &feed(),
            Temperature::new(300.0),
            controls(),
            || false,
        ),
        Err(GasEquilibriumChamberError::Invalid(_))
    ));
}

#[test]
fn structural_zeros_drop_only_unavailable_elements_and_retain_product_sources() {
    let (bank, chemistry) = fixture(
        vec![species("A", 0.01, 4.0, 0.0), species("B", 0.02, 3.5, 0.0)],
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2.0e-5, 0.15)],
        &["H", "O"],
        vec![vec![1, 0], vec![0, 1]],
    );
    let result = bank
        .equilibrium_chamber_hp(
            &chemistry,
            &point(300.0, 1e5),
            &feed(),
            Temperature::new(300.0),
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(result.equilibrium().species()[1].log_mole_fraction, None);
    assert_eq!(
        result.parcel().model().thermodynamics().components(),
        &chemistry.models()[..1]
    );
    assert_eq!(
        result.parcel().model().transport_materials(),
        &bank.transport_materials()[..1]
    );
    assert_eq!(
        result
            .parcel()
            .model()
            .thermodynamics()
            .mole_composition()
            .fractions(),
        &[1.0]
    );
    assert!(result.balance().enthalpy_residual_j.abs() < 1e-6);
}

#[test]
fn finite_logarithmic_underflow_is_refused_instead_of_dropping_reacting_species() {
    let a = species("A", 0.01, 4.0, 0.0);
    let b = species("B", 0.01, 4.0, 0.0);
    let mut card = b.regions()[0].card().clone();
    card.parameters.get_mut("a8").unwrap().value = 1000.0;
    let b = Nasa9StandardStateModelV1::new(
        b.species().clone(),
        b.molar_mass(),
        b.convention().clone(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap();
    let (bank, chemistry) = fixture(
        vec![a, b],
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2.0e-5, 0.15)],
        &["H"],
        vec![vec![1, 1]],
    );
    let result = bank.equilibrium_chamber_hp(
        &chemistry,
        &point(1000.0, 1e5),
        &feed(),
        Temperature::new(1000.0),
        controls(),
        || false,
    );
    assert_eq!(
        result,
        Err(GasEquilibriumChamberError::UnrepresentableTrace(id("A")))
    );
}

#[test]
fn final_transport_domain_refusal_preserves_inputs_after_successful_chemistry() {
    let (bank, chemistry) = isomers();
    let limited = fixture(
        chemistry.models().to_vec(),
        vec![
            transport_in_domain("A", 1.7e-5, 0.025, [200.0, 350.0], [1000.0, 1e7]),
            transport_in_domain("B", 2.0e-5, 0.15, [200.0, 350.0], [1000.0, 1e7]),
        ],
        &["H"],
        vec![vec![1, 1]],
    )
    .0;
    let old = limited.clone();
    assert!(
        bank.equilibrium_chamber_hp(
            &chemistry,
            &point(300.0, 3e6),
            &feed(),
            Temperature::new(300.0),
            controls(),
            || false,
        )
        .is_ok()
    );
    assert!(matches!(
        limited.equilibrium_chamber_hp(
            &chemistry,
            &point(300.0, 3e6),
            &feed(),
            Temperature::new(300.0),
            controls(),
            || false,
        ),
        Err(GasEquilibriumChamberError::State(_))
    ));
    assert_eq!(limited, old);
}

#[test]
fn cancellation_at_entry_during_chemistry_and_before_publication_is_atomic() {
    let (bank, chemistry) = isomers();
    let old = (bank.clone(), chemistry.clone());
    let mut total = 0;
    bank.equilibrium_chamber_hp(
        &chemistry,
        &point(300.0, 3e6),
        &feed(),
        Temperature::new(300.0),
        controls(),
        || {
            total += 1;
            false
        },
    )
    .unwrap();
    for stop in [1, total / 2, total] {
        let mut calls = 0;
        assert_eq!(
            bank.equilibrium_chamber_hp(
                &chemistry,
                &point(300.0, 3e6),
                &feed(),
                Temperature::new(300.0),
                controls(),
                || {
                    calls += 1;
                    calls == stop
                },
            ),
            Err(GasEquilibriumChamberError::Cancelled)
        );
        assert_eq!((&bank, &chemistry), (&old.0, &old.1));
    }
    for mass_tolerance_kg in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            bank.equilibrium_chamber_hp(
                &chemistry,
                &point(300.0, 3e6),
                &feed(),
                Temperature::new(300.0),
                GasEquilibriumChamberSolve {
                    mass_tolerance_kg,
                    ..controls()
                },
                || false,
            ),
            Err(GasEquilibriumChamberError::Invalid(_))
        ));
    }
}
